#include "turbulence/Turbulence.h"

double SA::_chi(double rho, double mu, double miubl){return rho*(miubl > 0.0 ? miubl : 0.0)/mu;}
