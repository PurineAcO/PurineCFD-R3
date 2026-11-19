#include "data_structure/Vec3.h"

namespace cc {

vec3 unit(const vec3& a) {return (1.0/a.norm())*a;}

}
