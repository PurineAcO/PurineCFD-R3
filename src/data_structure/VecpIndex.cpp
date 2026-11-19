#include "data_structure/Vecp.h"

namespace cc {

double& vecp::operator[](int i) {return (&rho)[i];}

double vecp::operator[](int i) const {return (&rho)[i];}

}
