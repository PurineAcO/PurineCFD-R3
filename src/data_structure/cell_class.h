#pragma once
#include "classconfig.h"

namespace cc {
// 格子
struct cell_class {
    int index = 0;                             // 编号
    int fcnt = 6;                              // 面数
    bool fnorm[6] = {};                        // 邻接面外法向标记
    cell_class* nei[6] = {};                   // 邻接网格指针
    face_class* faces[6] = {};                 // 邻接面指针
    double vol = 0.0;                          // 体积
    vec3 center;                               // 中心坐标

    vecp phy;                                   // 物理量
    struct otphy* otphy = nullptr;               // 引申物理量
    vec5 conser;                                // 守恒量
    vec5 lastconser;                            // 上一时间步守恒量
    vec5 convect;                               // 无粘对流项
    vec5 visflux;                               // 黏性通量
    vecgrad phgrad;                             // 物理量梯度

    // #if JSTchosen
    //     JST::dissipation diss;     // 耗散项
    // #elif ROEchosen
    //     ROE::dissipation diss;     // 耗散项
    // #endif

    #if SAchosen
        SA::turbulence* tur = NULL;   // 湍流
    #elif SSTchosen
        SST::turbulence* tur = NULL;  // 湍流
    #endif 

    #if LSCBchosen
        struct LSCBmatrix* LSCB = NULL; // LSCB梯度预处理矩阵
    #endif
    mat5 diag;            // 隐式对角块
    mat5tri conj;         // 无粘通量Jacobi
    double lambda = 0.0;  // 伪时间谱半径

    // 由守恒量还原基本物理量
    void prim();
    // 由基本物理量形成守恒量
    void form_conservative();

    cell_class() = default;
    cell_class(int index_, int f1_, int f2_, int f3_, int f4_, int f5_, int f6_,int fcnt_);
};

}