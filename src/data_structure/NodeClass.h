#pragma once
#include <cstddef>

#include "data_structure/Vec3.h"

namespace cc {

struct node_class {
    int index = 0;
    vec3 p;
    node_class() = default;
    node_class(int number_,double x_,double y_,double z_):index(number_),p(vec3(x_,y_,z_)){}
};

}
