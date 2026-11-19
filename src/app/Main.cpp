#include "app/App.h"
#include "config/Config.h"
#include "data_structure/DataStructure.h"
#include "explicit/Explicit.h"
#include "geometry/Geometry.h"
#include "implicit/Implicit.h"
#include "initialization/Initialization.h"
#include "io/Io.h"
#include "mesh/Mesh.h"
#include "physics/Physics.h"
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef PURINE_DATA_DIR
#define PURINE_DATA_DIR "run"
#endif

int main(){
    setvbuf(stdout, nullptr, _IONBF, 0);
#ifdef _OPENMP
    if(cc::threads > 0) {omp_set_num_threads(cc::threads);}
#endif
    const std::string meshpath = solver::mesh;
    const int maxstep = solver::max_step;
    const double cflmax = solver::cfl;
    const int turboff = solver::turboff;
    const std::string meshtag = std::filesystem::path(meshpath).parent_path().filename().string();
    const std::string autotag = std::string(solver::tag).empty()
        ? (meshtag.empty() ? std::string("case") : meshtag) +
          (dualtime::dt > 0.0 ? "-urans" : "-steady")
        : std::string(solver::tag);
    const std::string outdir = std::string(solver::outdir).empty()
        ? next_case_dir(PURINE_DATA_DIR, autotag)
        : std::string(solver::outdir);
    std::filesystem::create_directories(outdir);
    std::printf("output: %s\n", outdir.c_str());
    std::printf("solver: cfl %.2f  sweeps %d  conv_step %d  turboff %d\n",
                cflmax, solver::sweeps, solver::conv_step, turboff);
#ifdef _OPENMP
    std::printf("threads: %d (config %d, 0=auto)\n", omp_get_max_threads(), cc::threads);
#endif
    if(dualtime::dt > 0.0){
        std::printf("dualtime: dt %.4e  inner %d  tol %.1e  steady_iters %d  steady_cfl %.1f  seed %.4f\n",
                    dualtime::dt, dualtime::inner, dualtime::inner_tol, dualtime::steady_iters,
                    dualtime::steady_cfl, dualtime::seed);
    }

    if(!readmesh(meshpath.c_str())) {return 1;}
    allocate_fields();
    geometry_main();
#pragma omp parallel for schedule(static)
    allcell icell(i).tur->sad = distance_to_wall(icell(i));
    standard_initialization();
#pragma omp parallel for schedule(static)
    allcell icell(i).form_conservative();

    const cc::vecp farphy = from_boundary_to_phy();
    const double pref = farphy.rho*ph::R*farphy.T;
    std::printf("mesh: %d nodes, %d faces, %d cells\n", cc::node_num, cc::face_num, cc::cell_num);

    int outward = 0, inward = 0;
    for(const cc::bound& b : cc::bounds){
        if(b.type == cc::INTER) {continue;}
        for(const cc::face_class* f : b.fs){
            const cc::cell_class* c = cc::boundary_findcell(const_cast<cc::face_class*>(f));
            if(cc::dot(f->nor, f->mid - c->center) > 0.0) {outward++;} else {inward++;}
        }
    }
    std::printf("boundary normal: outward=%d inward=%d\n", outward, inward);


    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]{
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    FILE* hist = nullptr;
    const std::string histpath = outdir + "/history.csv";
    const std::string fieldpath = outdir + "/field.csv";
    const std::string wallpath = outdir + "/wall_cp.csv";
    const errno_t herr = fopen_s(&hist, histpath.c_str(), "w");
    if(hist == nullptr) {
        std::fprintf(stderr, "无法创建 %s (errno %d)\n", histpath.c_str(), (int)herr);
        return 1;
    }
    fprintf(hist, "phase,step,time,inner,res0,res,cl,cd\n");
    FILE* cpmid = nullptr;
    const errno_t cerr = fopen_s(&cpmid, (outdir + "/wall_cp_mid.csv").c_str(), "w");
    if(cpmid == nullptr) {
        std::fprintf(stderr, "无法创建 wall_cp_mid.csv (errno %d)\n", (int)cerr);
        fclose(hist);
        return 1;
    }
    fprintf(cpmid, "step,theta,cp\n");

    if(dualtime::seed > 0.0) {
        add_seed(dualtime::seed*std::hypot(farphy.u, farphy.v));
        std::printf("seed: %.4e (%.2f%% Uinf)\n", dualtime::seed, 100.0*dualtime::seed);
    }

    double cl = 0.0, cd = 0.0;
    int step = 0;
    if(dualtime::dt > 0.0){
        for(int it = 1; it <= dualtime::steady_iters;){
            const double cfl_nom = std::min(dualtime::steady_cfl, 2.0 + (dualtime::steady_cfl - 2.0)*it/200.0);
            save_state();
            const double r = solver::steady_scheme == 0
                ? impl::flow_step(0.0, cfl_nom*cfl_scale)
                : impl::flow_step_explicit(0.0);
            if(!turboff) {
                if(solver::steady_scheme == 0) {SA::turbulence_step(0.0, cfl_nom*cfl_scale);}
                else {SA::turbulence_step_explicit(0.0);}
            }
            const int nb = ph::bad_cells();
            if(nb > 0){
                printf("steady %6d 发散: %d 个坏单元, 回退上一流场, CFL x%.4f -> x%.4f\n",
                       it, nb, cfl_scale, cfl_scale*0.5);
                fflush(stdout);
                restore_state();
                cfl_scale *= 0.5;
                if(cfl_scale < 1e-3){printf("CFL 已降到下限仍不收敛, 放弃\n"); return 1;}
                continue;
            }
            if(cfl_scale < 1.0) {cfl_scale = std::min(1.0, cfl_scale*1.05);}
#pragma omp parallel for schedule(static)
            allcell icell(i).prim();
            if(it % solver::conv_step == 0 || it == dualtime::steady_iters){
                ph::wall_forces(cl, cd, farphy);
                fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", it, r, r, cl, cd);
                fflush(hist);
                printf("steady %6d  cfl %6.1f  res %.4e  Cl %+.5f  Cd %+.5f  clip %5d  [%.0fs]\n",
                       it, cfl_nom*cfl_scale, r, cl, cd, impl::n_clip, elapsed());
                fflush(stdout);
            }
            ++it;
        }
        if(dualtime::steady_iters > 0){
#pragma omp parallel for schedule(static)
            allcell icell(i).prim();
            if(!dump_field((outdir + "/field_steady.csv").c_str())) {printf("dump failed\n");}
        }
        double t = 0.0;
        for(step = 1; step <= maxstep; step++){
#pragma omp parallel for schedule(static)
            for(int i = 0; i < cc::cell_num; i++){
                cc::CellList[i].lastconser = cc::CellList[i].conser;
                cc::CellList[i].tur->miubl_former = cc::CellList[i].tur->miubl;
            }
            double r0 = 0.0, r = 0.0;
            int k = 0;
            const double cfl_cap = std::min(cflmax, 10.0);
            const double cfl_urans = std::min(cfl_cap, 2.0 + (cfl_cap - 2.0)*(step - 1)/200.0);
            while(k < dualtime::inner){
                r = solver::inner_scheme == 0
                    ? impl::flow_step(dualtime::dt, cfl_urans)
                    : impl::flow_step_explicit(dualtime::dt);
                if(!turboff) {
                    if(solver::inner_scheme == 0) {SA::turbulence_step(dualtime::dt, cfl_urans);}
                    else {SA::turbulence_step_explicit(dualtime::dt);}
                }
#pragma omp parallel for schedule(static)
                allcell icell(i).prim();
                if(++k == 1) {r0 = r;} else if(r < dualtime::inner_tol*r0) {break;}
            }
            t += dualtime::dt;
            ph::wall_forces(cl, cd, farphy);
            write_wall_cp_mid(cpmid, step, pref, farphy);
            fprintf(hist, "U,%d,%.8e,%d,%.6e,%.6e,%.8f,%.8f\n", step, t, k, r0, r, cl, cd);
            fflush(hist);
            if(!std::isfinite(cl) || !std::isfinite(cd) || ph::bad_cells() > 0){
                printf("step %d t=%.6e inner=%d res0=%.6e res=%.6e diverge\n", step, t, k, r0, r);
                break;
            }
            if(step % solver::conv_step == 0 || step == 1){
                double numax = 0.0;
                for(int i = 0; i < cc::cell_num; i++) {numax = std::max(numax, cc::CellList[i].tur->miubl);}
                printf("step %6d  t %.5e  inner %2d  res %.2e -> %.2e  Cl %+.5f  Cd %+.5f  clip %5d  miubl_max %.3e  [%.0fs]\n",
                       step, t, k, r0, r, cl, cd, impl::n_clip, numax, elapsed());
                fflush(stdout);
            }
            if(step % solver::field_interval == 0 && step >= solver::field_start_step){
#pragma omp parallel for schedule(static)
                allcell icell(i).prim();
                char buf[512];
                std::snprintf(buf, sizeof(buf), "%s/field_%05d.csv", outdir.c_str(), step);
                dump_field(buf);
            }
        }
    }else{
        for(int step_it = 1; step_it <= maxstep;){
            const double cfl_nom = std::min(cflmax, 2.0 + (cflmax - 2.0)*step_it/200.0);
            save_state();
            const double res_now = solver::steady_scheme == 0
                ? impl::flow_step(0.0, cfl_nom*cfl_scale)
                : impl::flow_step_explicit(0.0);
            if(!turboff) {
                if(solver::steady_scheme == 0) {SA::turbulence_step(0.0, cfl_nom*cfl_scale);}
                else {SA::turbulence_step_explicit(0.0);}
            }
            const int nb = ph::bad_cells();
            if(nb > 0){
                printf("step %6d 发散: %d 个坏单元, 回退上一流场, CFL x%.4f -> x%.4f\n",
                       step_it, nb, cfl_scale, cfl_scale*0.5);
                fflush(stdout);
                restore_state();
                cfl_scale *= 0.5;
                if(cfl_scale < 1e-3){printf("CFL 已降到下限仍不收敛, 放弃\n"); break;}
                continue;
            }
            if(cfl_scale < 1.0) {cfl_scale = std::min(1.0, cfl_scale*1.05);}
            step = step_it;
#pragma omp parallel for schedule(static)
            allcell icell(i).prim();
            ph::wall_forces(cl, cd, farphy);
            fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", step, res_now, res_now, cl, cd);
            fflush(hist);
            write_wall_cp_mid(cpmid, step, pref, farphy);
            if(step % solver::conv_step == 0){
                double numax = 0.0;
                for(int i = 0; i < cc::cell_num; i++) {numax = std::max(numax, cc::CellList[i].tur->miubl);}
                printf("step %6d  cfl %6.1f  res %.4e  Cl %+.5f  Cd %+.5f  clip %5d  miubl_max %.4e  [%.0fs]\n",
                       step, cfl_nom*cfl_scale, res_now, cl, cd, impl::n_clip, numax, elapsed());
                fflush(stdout);
            }
            ++step_it;
        }
    }
    if(step > maxstep) {step = maxstep;}
#pragma omp parallel for schedule(static)
    allcell icell(i).prim();
    ph::wall_forces(cl, cd, farphy);
    fclose(hist);
    if(cpmid != nullptr) {fclose(cpmid);}
    dump_field(fieldpath.c_str());
    write_wall_cp(wallpath.c_str(), pref, farphy);
    printf("step %d  Cl=%.6f Cd=%.6f clip=%d  [%.0fs]\n", step, cl, cd, impl::n_clip, elapsed());
    return 0;
}
