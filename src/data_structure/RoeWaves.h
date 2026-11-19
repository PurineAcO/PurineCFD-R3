#pragma once

#include "data_structure/Vec3.h"
#include "data_structure/Vec5.h"

namespace ROE {
struct dissipation;
}

struct roe_waves {
    cc::vec3 n;
    double alpha[5];
    double lambda[5];
    cc::vec5 r[5];
};
