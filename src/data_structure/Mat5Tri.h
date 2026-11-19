#pragma once

#include "data_structure/Mat5.h"

namespace cc {

struct mat5tri {
    mat5 x, y, z;
    mat5tri() = default;
    mat5tri(mat5 x_, mat5 y_, mat5 z_) : x(x_), y(y_), z(z_) {}
};

}
