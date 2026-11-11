#pragma once
#include "classconfig.h"
#include <cstddef>


namespace cc {

// 面隐式Jacobian
struct implict {
    mat5 Bf;  // 冻结梯度粘性Jacobi
    mat5 Tf;  // 梯度粘性Jacobi
    mat5 Ab;  // Roe通量Jacobi
};

// 面
struct face_class {
    int index = 0;                // 面编号
    short type = INTER;           // 面类型
    short pcnt = 4;               // 节点数
    node_class* node[4] = {};     // 面邻接点指针(按右手定则排序)
    vec3 mid = {};                // 面中点坐标
    vec3 nor = {};                // 单位法向量
    bool outer = false;           // 面法向与nei[0]->nei[1]同向时为true
    double area = 0.0;            // 面面积
    cell_class* nei[2] = {};      // 面邻接网格指针
    vecp phy;                     // 物理量
    vecgrad phgrad;               // 面法向梯度
    struct otphy* otphy = nullptr; // 引申物理量
    double un = 0.0;              // 面法向速度
    #if SAchosen
        SA::turbulence* tur = NULL;   // 湍流
    #elif SSTchosen
        SST::turbulence* tur = NULL;  // 湍流
    #endif
    vecp phynei[2];               // 左右插值物理量
    struct otphy* otnei[2] = {};   // 左右引申物理量
    implict imats;                // 非对角矩阵
    vec5 convect;                 // 无粘通量
    vec5 visflux;                 // 黏性通量

    // 标量投影到面法向
    double toface_jacobi(double F, double G, double K) const;
    // 通量投影到面法向
    vec5 toface_jacobi(const vec5& F, const vec5& G, const vec5& K) const;
    // 矩阵投影到面法向
    mat5 toface_jacobi(const mat5tri& m) const;

    face_class() = default;
    face_class(int index_, int p1_, int p2_, int p3_, int p4_, int c1_, int c2_, short type_,short pcnt_);
};

}