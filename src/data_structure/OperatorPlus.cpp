#include "data_structure/Mat5.h"
#include "data_structure/Vec3.h"
#include "data_structure/Vec5.h"
#include "data_structure/Vecgrad.h"
#include "data_structure/Vecp.h"

namespace cc {

vec3 operator+(vec3 a, const vec3& b) {return a += b;}

vecp operator+(vecp a, const vecp& b) {return a += b;}

vec5 operator+(vec5 a, const vec5& b) {return a += b;}

vecgrad operator+(const vecgrad& a, const vecgrad& b) {
    return vecgrad(a.rhograd + b.rhograd, a.ugrad + b.ugrad, a.vgrad + b.vgrad,
                   a.wgrad + b.wgrad, a.Tgrad + b.Tgrad);
}

mat5 operator+(mat5 m, const mat5& o) { return m += o; }

}
