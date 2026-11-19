#include "turbulence/Turbulence.h"

double SA::_g(double r){return r + SA::Cw2*(r*r*r*r*r*r - r);}
