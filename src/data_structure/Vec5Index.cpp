#include "data_structure/Vec5.h"

namespace cc {

double& vec5::operator[](int i) {return (&c)[i];}

double vec5::operator[](int i) const {return (&c)[i];}

}
