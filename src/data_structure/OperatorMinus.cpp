#include "data_structure/Mat5.h"
#include "data_structure/Vec3.h"
#include "data_structure/Vec5.h"
#include "data_structure/Vecp.h"

namespace cc {

vec3 operator-(vec3 a, const vec3& b) {return a -= b;}

vecp operator-(vecp a, const vecp& b) {return a -= b;}

vec5 operator-(vec5 a, const vec5& b) {return a -= b;}

mat5 operator-(mat5 m, const mat5& o) {
    for(int r = 0; r < 5; r++)
        for(int c = 0; c < 5; c++) m.a[r][c] -= o.a[r][c];
    return m;
}

}
