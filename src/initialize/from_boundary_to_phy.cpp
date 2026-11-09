#include "initialization.h"
#include "../physic/physic.h"
#include <cmath>

cc::vecp from_boundary_to_phy(){
    const double a = std::sqrt(ph::gamma*ph::R*cc::FAR_config.T);
    const double v = cc::FAR_config.Ma*a;
    return cc::vecp(cc::FAR_config.p/(ph::R*cc::FAR_config.T),
                    v*std::cos(cc::FAR_config.AOA),v*std::sin(cc::FAR_config.AOA),0.0,
                    cc::FAR_config.T);
}