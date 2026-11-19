#pragma once
#include <string>
#include <vector>

#include "data_structure/BoundaryTypes.h"

namespace cc {

struct face_class;

struct bound {
    std::string name;
    short type = INTER;
    std::vector<face_class*> fs;
    bound() = default;
    bound(std::string name_,short type_):name(name_),type(type_) {}
};

}
