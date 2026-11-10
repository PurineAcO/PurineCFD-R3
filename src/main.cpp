#include "data_structure/cell_class.h"
#include "data_structure/classconfig.h"
#include "data_structure/face_class.h"
#include "geometry/geometry.h"
#include "implict.hpp"
#include "initialize/initialization.h"
#include "io.hpp"
#include "mesh/readmesh.h"
#include "physic/physic.h"
#include "residual.hpp"
#include "timarch.hpp"
#include "udf/udf.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#ifndef PURINE_DATA_DIR
#define PURINE_DATA_DIR "run"
#endif

namespace {

struct urans_cfg {
    double dt = 1.5e-3;
    int inner = 20;
    double inner_tol = 1e-3;
    int steady_iters = 200;
    double steady_cfl = 20.0;
    double seed = 0.02;
};

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

int main(int argc, char** argv){
    setvbuf(stdout, nullptr, _IONBF, 0);
    const char* meshpath = argc > 1 ? argv[1] : "cases/cyc/mesh.txt";
    const int maxstep = argc > 2 ? std::atoi(argv[2]) : 2000;
    const double cflmax = argc > 3 ? std::atof(argv[3]) : 1.0;
    const char mode = argc > 4 ? argv[4][0] : 'i';
    const int turboff = argc > 5 ? std::atoi(argv[5]) : 0;
    if(argc > 6) {impl::sweeps = std::atoi(argv[6]);}
    if(argc > 7) {res::conv_step = std::atoi(argv[7]);}
    urans_cfg ur;
    if(argc > 9) {ur.dt = std::atof(argv[9]);}
    if(argc > 10) {ur.inner = std::atoi(argv[10]);}
    if(argc > 11) {ur.steady_iters = std::atoi(argv[11]);}
    if(argc > 12) {ur.steady_cfl = std::atof(argv[12]);}
    if(argc > 13) {ur.seed = std::atof(argv[13]);}
    const std::string meshtag = std::filesystem::path(meshpath).parent_path().filename().string();
    const std::string tag = argc > 14 ? argv[14]
        : (meshtag.empty() ? std::string("case") : meshtag) + (ur.dt > 0.0 ? "-urans" : "-steady");
    const std::string outdir = argc > 8 ? argv[8] : next_case_dir(PURINE_DATA_DIR, tag);
    std::filesystem::create_directories(outdir);
    std::printf("output: %s\n", outdir.c_str());

    if(!readmesh(meshpath)) {return 1;}
    allocate_fields();
    geometry_main();
    allcell icell(i).tur->sad = distance_to_wall(icell(i));
    standard_initialization();
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

    if(ur.seed > 0.0) {
        add_seed(ur.seed*std::hypot(farphy.u, farphy.v));
        std::printf("seed: %.4e (%.2f%% Uinf)\n", ur.seed, 100.0*ur.seed);
    }

    double cl = 0.0, cd = 0.0;
    int step = 0;
    if(ur.dt > 0.0){
        for(int it = 1; it <= ur.steady_iters; it++){
            const double cfl = std::min(ur.steady_cfl, 2.0 + (ur.steady_cfl - 2.0)*it/200.0);
            const double r = impl::flow_step(0.0, cfl);
            if(!turboff) {impl::sa_step(0.0, cfl);}
            allcell icell(i).prim();
            if(it % res::conv_step == 0 || it == ur.steady_iters){
                ph::wall_forces(cl, cd, farphy);
                fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", it, r, r, cl, cd);
                fflush(hist);
                printf("steady %6d  cfl %6.1f  res %.4e  Cl %+.5f  Cd %+.5f  clip %5d  [%.0fs]\n",
                       it, cfl, r, cl, cd, impl::n_clip, elapsed());
                fflush(stdout);
            }
        }
        if(ur.steady_iters > 0){
            allcell icell(i).prim();
            if(!dump_field((outdir + "/field_steady.csv").c_str())) {printf("dump failed\n");}
        }
        double t = 0.0;
        for(step = 1; step <= maxstep; step++){
            for(int i = 0; i < cc::cell_num; i++){
                cc::CellList[i].lastconser = cc::CellList[i].conser;
                cc::CellList[i].tur->miubl_former = cc::CellList[i].tur->miubl;
            }
            double r0 = 0.0, r = 0.0;
            int k = 0;
            while(k < ur.inner){
                r = impl::flow_step(ur.dt, cflmax);
                if(!turboff) {impl::sa_step(ur.dt, cflmax);}
                allcell icell(i).prim();
                if(++k == 1) {r0 = r;} else if(r < ur.inner_tol*r0) {break;}
            }
            t += ur.dt;
            ph::wall_forces(cl, cd, farphy);
            write_wall_cp_mid(cpmid, step, pref, farphy);
            fprintf(hist, "U,%d,%.8e,%d,%.6e,%.6e,%.8f,%.8f\n", step, t, k, r0, r, cl, cd);
            fflush(hist);
            if(!std::isfinite(cl) || !std::isfinite(cd) || ph::bad_cells() > 0){
                printf("step %d t=%.6e inner=%d res0=%.6e res=%.6e diverge\n", step, t, k, r0, r);
                break;
            }
            if(step % res::conv_step == 0 || step == 1){
                double numax = 0.0;
                for(int i = 0; i < cc::cell_num; i++) {numax = std::max(numax, cc::CellList[i].tur->miubl);}
                printf("step %6d  t %.5e  inner %2d  res %.2e -> %.2e  Cl %+.5f  Cd %+.5f  clip %5d  miubl_max %.3e  [%.0fs]\n",
                       step, t, k, r0, r, cl, cd, impl::n_clip, numax, elapsed());
                fflush(stdout);
            }
            if(step % 100 == 0){
                allcell icell(i).prim();
                char buf[512];
                std::snprintf(buf, sizeof(buf), "%s/field_%05d.csv", outdir.c_str(), step);
                dump_field(buf);
            }
        }
    }else{
        for(step = 1; step <= maxstep; step++){
            const double cfl = std::min(cflmax, 2.0 + (cflmax - 2.0)*step/200.0);
            double res_now = 0.0;
            if(mode == 'i'){
                res_now = impl::flow_step(0.0, cfl);
                if(!turboff) {impl::sa_step(0.0, cfl);}
                allcell icell(i).prim();
            }else{
                fatime::CFL = cfl;
                one_rans(0.0, false);
                allcell icell(i).prim();
                calculate_residual();
                res_now = residual_drop();
            }
            ph::wall_forces(cl, cd, farphy);
            fprintf(hist, "S,%d,0,1,%.6e,%.6e,%.8f,%.8f\n", step, res_now, res_now, cl, cd);
            fflush(hist);
            write_wall_cp_mid(cpmid, step, pref, farphy);
            if(step % res::conv_step == 0){
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
    allcell icell(i).prim();
    ph::wall_forces(cl, cd, farphy);
    fclose(hist);
    if(cpmid != nullptr) {fclose(cpmid);}
    dump_field(fieldpath.c_str());
    write_wall_cp(wallpath.c_str(), pref, farphy);
    printf("step %d  Cl=%.6f Cd=%.6f clip=%d  [%.0fs]\n", step, cl, cd, impl::n_clip, elapsed());
    return 0;
}
