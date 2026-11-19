#pragma once

#include "data_structure/Vec3.h"

namespace cc {

struct LSCBmatrix {
    double Axx = 0.0;
    double Axy = 0.0;
    double Axz = 0.0;
    double Ayy = 0.0;
    double Ayz = 0.0;
    double Azz = 0.0;
    double wi[6] = {};
    double invA[3][3] = {};
    LSCBmatrix() = default;
    // 清空
    void clear();
    // 形成几何逆
    void form_invA();
    // 解梯度方程
    cc::vec3 solve(cc::vec3 dphi);
};

inline void LSCBmatrix::clear() {
    Axx = Ayy = Azz = Axy = Axz = Ayz = 0.0;
}

inline void LSCBmatrix::form_invA() {
    invA[0][0] = Ayy*Azz - Ayz*Ayz;
    invA[1][1] = Axx*Azz - Axz*Axz;
    invA[2][2] = Axx*Ayy - Axy*Axy;
    invA[0][1] = invA[1][0] = Axz*Ayz - Axy*Azz;
    invA[0][2] = invA[2][0] = Axy*Ayz - Axz*Ayy;
    invA[1][2] = invA[2][1] = Axy*Axz - Ayz*Axx;
    const double det = Axx*invA[0][0] + Axy*invA[0][1] + Axz*invA[0][2];
    for(int r=0;r<3;r++)
        for(int s=0;s<3;s++) invA[r][s] /= det;
}

inline cc::vec3 LSCBmatrix::solve(cc::vec3 dphi) {
    return cc::vec3(invA[0][0]*dphi.x + invA[0][1]*dphi.y + invA[0][2]*dphi.z,
                    invA[1][0]*dphi.x + invA[1][1]*dphi.y + invA[1][2]*dphi.z,
                    invA[2][0]*dphi.x + invA[2][1]*dphi.y + invA[2][2]*dphi.z);
}

}
