#include "turbulence/Turbulence.h"

double SA::_ft2(double chi){return SA::Ct3*std::exp(-SA::Ct4*chi*chi);}
