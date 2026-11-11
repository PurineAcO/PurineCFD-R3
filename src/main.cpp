#include "data_structure/cell_class.h"
#include "data_structure/classconfig.h"
#include "data_structure/face_class.h"
#include "physic/physic.h"
#include "udf/udf.h"
#include "DensityBased/implict.h"
#include "geometry/geometry.h"
#include "initialize/initialization.h"
#include "io.hpp"
#include "mesh/readmesh.h"
#include <algorithm>
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

namespace {

std::string next_case_dir(const std::string& root, const std::string& tag){
    int next = 1;
    std::error_code ec;
    if(std::filesystem::is_directory(root, ec)){
        for(const std::filesystem::directory_entry& e : std::filesystem::directory_iterator(root, ec)){
            const std::string name = e.path().filename().string();
            if(name.compare(0, 3, "3D-") != 0) {continue;}
            const std::size_t dash = name.find('-', 3);
            if(dash == std::string::npos) {continue;}
            const int v = std::atoi(name.substr(3, dash - 3).c_str());
            next = std::max(next, v + 1);
        }
    }
    return root + "/3D-" + std::to_string(next) + "-" + tag;
}

void add_seed(double amp){
    for(cc::cell_class& c : cc::CellList){
        const double r = std::hypot(c.center.x, c.center.y);
        if(r < 0.5 || r > 3.0) {continue;}
        const double wgt = amp*std::exp(-r*r/4.0);
        c.phy.u -= wgt*c.center.y;
        c.phy.v += wgt*c.center.x;
        *c.otphy = ph::form_otphy(c.phy);
        c.form_conservative();
    }
}

}

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
    fopen_s(&hist, histpath.c_str(), "w");
    fprintf(hist, "phase,step,time,inner,res0,res,cl,cd\n");
    FILE* cpmid = nullptr;
    fopen_s(&cpmid, (outdir + "/wall_cp_mid.csv").c_str(), "w");
    fprintf(cpmid, "step,theta,cp\n");

    if(dualtime::seed > 0.0) {
        add_seed(dualtime::seed*std::hypot(farphy.u, farphy.v));
        std::printf("seed: %.4e (%.2f%% Uinf)\n", dualtime::seed, 100.0*dualtime::seed);
    }

    double cl = 0.0, cd = 0.0;
    int step = 0;
    if(dualtime::dt > 0.0){
        for(int it = 1; it <= dualtime::steady_iters; it++){
            const double cfl = std::min(dualtime::steady_cfl, 2.0 + (dualtime::steady_cfl - 2.0)*it/200.0);
            const double r = impl::flow_step(0.0, cfl);
            if(!turboff) {SA::turbulence_step(0.0, cfl);}
#pragma omp parallel for schedule(static)
            allcell icell(i).prim();
            if(it % solver::conv_step == 0 || it == dualtime::steady_iters){
                ph::wall_forces(cl, cd, farphy);
                fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", it, r, r, cl, cd);
                fflush(hist);
                printf("steady %6d  cfl %6.1f  res %.4e  Cl %+.5f  Cd %+.5f  clip %5d  [%.0fs]\n",
                       it, cfl, r, cl, cd, impl::n_clip, elapsed());
                fflush(stdout);
            }
            if(ph::bad_cells() > 0) {printf("steady %d diverge\n", it); return 1;}
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
            while(k < dualtime::inner){
                r = impl::flow_step(dualtime::dt, cflmax);
                if(!turboff) {SA::turbulence_step(dualtime::dt, cflmax);}
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
            if(step % solver::field_interval == 0){
#pragma omp parallel for schedule(static)
                allcell icell(i).prim();
                char buf[512];
                std::snprintf(buf, sizeof(buf), "%s/field_%05d.csv", outdir.c_str(), step);
                dump_field(buf);
            }
        }
    }else{
        for(step = 1; step <= maxstep; step++){
            const double cfl = std::min(cflmax, 2.0 + (cflmax - 2.0)*step/200.0);
            const double res_now = impl::flow_step(0.0, cfl);
            if(!turboff) {SA::turbulence_step(0.0, cfl);}
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
                       step, cfl, res_now, cl, cd, impl::n_clip, numax, elapsed());
                fflush(stdout);
            }
            if(ph::bad_cells() > 0) {printf("step %d diverge\n", step); break;}
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
