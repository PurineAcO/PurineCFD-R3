#include "data_structure/Mat5.h"
#include "data_structure/Vec3.h"
#include "data_structure/Vec5.h"
#include "data_structure/Vecgrad.h"
#include "data_structure/Vecp.h"

namespace cc {

vec3 operator*(vec3 a, double s) {return a *= s;}

vec3 operator*(double s, vec3 a) {return a *= s;}

vecp operator*(vecp a, double s) {return a *= s;}

vecp operator*(double s, vecp a) {return a *= s;}

vec5 operator*(vec5 a, double s) {return a *= s;}

vec5 operator*(double s, vec5 a) {return a *= s;}

vecp operator*(const vecgrad& a, const vec3& r) {
    return vecp(dot(a.rhograd, r), dot(a.ugrad, r), dot(a.vgrad, r), dot(a.wgrad, r),
                dot(a.Tgrad, r));
}

vecgrad operator*(const vecp& phyf, const vec3& fnor) {
    return vecgrad(phyf.rho * fnor, phyf.u * fnor, phyf.v * fnor, phyf.w * fnor, phyf.T * fnor);
}

vecgrad operator*(double a, const vecgrad& grad) {
    return vecgrad(a * grad.rhograd, a * grad.ugrad, a * grad.vgrad, a * grad.wgrad,
                   a * grad.Tgrad);
}

mat5 operator*(mat5 m, double k) { return m *= k; }

mat5 operator*(double k, mat5 m) { return m *= k; }

vec5 operator*(const mat5& m, const vec5& v) {
    return vec5(m.a[0][0] * v.c + m.a[0][1] * v.x + m.a[0][2] * v.y + m.a[0][3] * v.z + m.a[0][4] * v.e,
                m.a[1][0] * v.c + m.a[1][1] * v.x + m.a[1][2] * v.y + m.a[1][3] * v.z + m.a[1][4] * v.e,
                m.a[2][0] * v.c + m.a[2][1] * v.x + m.a[2][2] * v.y + m.a[2][3] * v.z + m.a[2][4] * v.e,
                m.a[3][0] * v.c + m.a[3][1] * v.x + m.a[3][2] * v.y + m.a[3][3] * v.z + m.a[3][4] * v.e,
                m.a[4][0] * v.c + m.a[4][1] * v.x + m.a[4][2] * v.y + m.a[4][3] * v.z + m.a[4][4] * v.e);
}

}
