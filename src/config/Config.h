#pragma once

namespace boundary_settings {

struct FARconfig {
    double p;
    double Ma;
    double T;
    double AOA;
    double chi;
};

}

namespace cc {

#define chose2D 1
#define SSTchosen 0
#define SAchosen 1
#define ROEchosen 1
#define JSTchosen 0
#define LSCBchosen 1
inline int threads = 32;
#if chose2D
inline boundary_settings::FARconfig FAR_config = {63585.0302134262,0.72,300.0,4.5,3.0};
#else
inline boundary_settings::FARconfig FAR_config = {22.9,0.2,300,0.0,0.1};
#endif

}

namespace solver {
#if chose2D
inline const char* mesh = "cases/naca0012/mesh.txt";
#else
inline const char* mesh = "cases/cyc/mesh.txt";
#endif
#if chose2D
inline int max_step = 10000;
#else
inline int max_step = 100;
#endif
#if chose2D
inline double cfl = 100.0;
#else
inline double cfl = 50.0;
#endif
inline int turboff = 0;
inline int sweeps = 4;
inline int inner_scheme = 0;
inline int steady_scheme = 0;
inline double exp_cfl = 0.3;
inline int conv_step = 20;
inline int field_interval = 2;
inline int field_start_step = 100000000;
inline const char* outdir = "D:/PurineCFD-R3/run/NACA0012-Ma0.72-a4.5-Re1e7";
#if chose2D
inline const char* tag = "naca0012";
#else
inline const char* tag = "cyc-urans-rf";
#endif
}

namespace dualtime {
inline double dt = 0.0;
inline int inner = 40;
inline double inner_tol = 1e-3;
inline int steady_iters = 0;
inline double steady_cfl = 100.0;
inline double seed = 0.0;
}
