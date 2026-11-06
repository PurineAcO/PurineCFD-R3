#include "classconfig.hpp"
#include "config.hpp"
#include "geometry.hpp"
#include "implict.hpp"
#include "initialize.hpp"
#include "io.hpp"
#include "physic.hpp"
#include "readmesh.hpp"
#include "residual.hpp"
#include "timarch.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv){
    const char* meshpath = argc > 1 ? argv[1] : "cases/cyc/mesh.txt";
    const int maxstep = argc > 2 ? std::atoi(argv[2]) : 2000;
    const double cflmax = argc > 3 ? std::atof(argv[3]) : 20.0;
    const int turboff = argc > 4 ? std::atoi(argv[4]) : 0;
    const char mode = argc > 5 ? argv[5][0] : 'e';
    if(!readmesh(meshpath)) {return 1;}
    geometry_main();
    standard_initialization();
    allcell icell(i).form_conservative();
    const cc::vecp farphy = from_boundary_to_phy();
    const double pref = farphy.rho*cc::R*farphy.T;

    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&]{
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };

    FILE* hist = nullptr;
    fopen_s(&hist, "cases/cyc/history.csv", "w");
    fprintf(hist, "step,res,cl,cd\n");

    double cl = 0.0, cd = 0.0;
    int step = 0;
    for(step = 0; step < maxstep; step++){
        const double cfl = std::min(cflmax, 2.0 + (cflmax - 2.0)*step/200.0);
        if(mode == 'i'){
            impl::flow_step(0.0, cfl);
            if(!turboff) {impl::turbulence_step(0.0, cfl);}
        }else{
            one_rans(0.0, false);
        }
        allcell icell(i).prim();
        calculate_residual();
        if(step % res::conv_step == 0){
            wall_forces(cl, cd, farphy);
            report_residual(step, cl, cd, elapsed());
            fprintf(hist, "%d,%.9e,%.9e,%.9e\n", step, residual_drop(), cl, cd);
            fflush(hist);
        }
        if(residual_drop() < res::drop_target) {break;}
        if(bad_cells() > 0) {printf("step %d 发散\n", step); break;}
    }
    allcell icell(i).prim();
    wall_forces(cl, cd, farphy);
    report_residual(step, cl, cd, elapsed());
    fprintf(hist, "%d,%.9e,%.9e,%.9e\n", step, residual_drop(), cl, cd);
    fclose(hist);
    dump_field("cases/cyc/field.csv");
    write_wall_cp("cases/cyc/wall_cp.csv", pref, farphy);
    printf("Cl=%.5f Cd=%.5f clip=%d\n", cl, cd, impl::n_clip);
    return 0;
}
