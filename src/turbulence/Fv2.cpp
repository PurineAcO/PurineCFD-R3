#include "turbulence/Turbulence.h"

double SA::_fv2(double chi){return 1.0 - chi/(1.0 + chi*_fv1(chi));}
