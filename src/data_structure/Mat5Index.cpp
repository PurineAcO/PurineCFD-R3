#include "data_structure/Mat5.h"

namespace cc {

double* mat5::operator[](int r) { return a[r]; }

const double* mat5::operator[](int r) const { return a[r]; }

}
