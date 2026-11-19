#pragma once

#include "data_structure/Vec3.h"

namespace SA {

struct turbulence {
    double miubl = 0.0;
    double miubl_former = 0.0;
    double sad = 0.0;
    cc::vec3 miublgrad;
    double mueff = 0.0;
    double chi = 0.0;
    double fv1 = 0.0;
    double ft2 = 0.0;
    double fv2 = 0.0;
    double Omega = 0.0;
    double Sbl = 0.0;
    double S2 = 0.0;
    double r = 0.0;
    double g = 0.0;
    double fw = 0.0;
    double diag = 0.0;
};

}

namespace SST {
struct turbulence;
}
