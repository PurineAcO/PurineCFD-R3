#include "initialization/Initialization.h"
#include "data_structure/DataStructure.h"
#include "physics/Physics.h"
#include <cmath>

cc::vecp from_boundary_to_phy(){
    const double a = std::sqrt(ph::gamma*ph::R*cc::FAR_config.T);
    const double v = cc::FAR_config.Ma*a;
    const double aoa = cc::FAR_config.AOA*(std::acos(-1.0)/180.0);
    return cc::vecp(cc::FAR_config.p/(ph::R*cc::FAR_config.T),
                    v*std::cos(aoa),v*std::sin(aoa),0.0,
                    cc::FAR_config.T);
}
