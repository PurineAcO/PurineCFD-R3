#include "turbulence/Turbulence.h"

double SA::_fv1(double chi){return chi*chi*chi/(chi*chi*chi + SA::Cv1*SA::Cv1*SA::Cv1);}
