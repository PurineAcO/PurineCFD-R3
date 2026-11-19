#pragma once

#include "data_structure/Vecp.h"

namespace cc {

struct otphy {
    double a = 0.0;
    double p = 0.0;
    double e = 0.0;
    double un = 0.0;
    double mu = 0.0;
    otphy() = default;
    otphy(double a_, double p_, double e_);
    otphy(const cc::vecp& phy);
};

}
