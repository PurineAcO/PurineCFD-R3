#pragma once
#include <cstddef>

#include "data_structure/BoundaryTypes.h"
#include "data_structure/ClassConfig.h"
#include "data_structure/Implict.h"
#include "data_structure/Mat5Tri.h"
#include "data_structure/Vec3.h"
#include "data_structure/Vec5.h"
#include "data_structure/Vecgrad.h"
#include "data_structure/Vecp.h"

namespace cc {
struct otphy;
}
namespace SA { struct turbulence; }
namespace SST { struct turbulence; }

namespace cc {

struct face_class {
    int index = 0;
    short type = INTER;
    short pcnt = 4;
    node_class* node[4] = {};
    vec3 mid = {};
    vec3 nor = {};
    bool outer = false;
    double area = 0.0;
    cell_class* nei[2] = {};
    vecp phy;
    vecgrad phgrad;
    struct otphy* otphy = nullptr;
    double un = 0.0;
    #if SAchosen
    SA::turbulence* tur = NULL;
    #elif SSTchosen
    SST::turbulence* tur = NULL;
    #endif
    vecp phynei[2];
    struct otphy* otnei[2] = {};
    implict imats;
    vec5 convect;
    vec5 visflux;

    // 标量投影到面法向
    double toface_jacobi(double F, double G, double K) const;
    // 通量投影到面法向
    vec5 toface_jacobi(const vec5& F, const vec5& G, const vec5& K) const;
    // 矩阵投影到面法向
    mat5 toface_jacobi(const mat5tri& m) const;

    face_class() = default;
    face_class(int index_, int p1_, int p2_, int p3_, int p4_, int c1_, int c2_,
               short type_, short pcnt_);
};

}
