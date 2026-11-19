#include "turbulence/Turbulence.h"

double SA::_fw(double g){
    constexpr double cw3_sixth = SA::Cw3*SA::Cw3*SA::Cw3*SA::Cw3*SA::Cw3*SA::Cw3;
    const double g_squared = g*g;
    const double g_sixth = g_squared*g_squared*g_squared;
    return g*std::pow((1.0 + cw3_sixth)/(cw3_sixth + g_sixth), 1.0/6.0);
}
