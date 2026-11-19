#pragma once
#include <cstddef>

#include "data_structure/ClassConfig.h"
#include "data_structure/Mat5.h"
#include "data_structure/Mat5Tri.h"
#include "data_structure/Vec3.h"
#include "data_structure/Vec5.h"
#include "data_structure/Vecgrad.h"
#include "data_structure/Vecp.h"

namespace cc {
struct otphy;
struct LSCBmatrix;
}
namespace SA { struct turbulence; }
namespace SST { struct turbulence; }

namespace cc {

struct cell_class {
    int index = 0;
    int fcnt = 6;
    bool fnorm[6] = {};
    cell_class* nei[6] = {};
    face_class* faces[6] = {};
    double vol = 0.0;
    vec3 center;

    vecp phy;
    struct otphy* otphy = nullptr;
    vec5 conser;
    vec5 lastconser;
    vec5 convect;
    vec5 visflux;
    vecgrad phgrad;

    #if SAchosen
    SA::turbulence* tur = NULL;
    #elif SSTchosen
    SST::turbulence* tur = NULL;
    #endif

    #if LSCBchosen
    struct LSCBmatrix* LSCB = NULL;
    #endif

    mat5 diag;
    mat5tri conj;
    double lambda = 0.0;

    // 由守恒量还原基本物理量
    void prim();
    // 由基本物理量形成守恒量
    void form_conservative();

    cell_class() = default;
    cell_class(int index_, int f1_, int f2_, int f3_, int f4_, int f5_, int f6_, int fcnt_);
};

}
