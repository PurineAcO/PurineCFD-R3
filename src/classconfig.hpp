#pragma once
#include <string>
#include <vector>

#include "config.hpp"

#define allface(cell) for (int i = 0; i < (cell).ecnt; i++)
#define allcell for (int i = 0; i < cc::cell_num; i++)
#define icell(i) cc::CellList[i]
#define allfac for (int i = 0; i < cc::face_num; i++)
#define iface(i) cc::FaceList[i]

// 头
namespace cc {

struct cell_class; // 网格
struct face_class; // 面
struct node_class; // 节点
struct bound {
    std::string name;             // 名称
    short type = INTER;           // 类型
    std::vector<face_class*> fs;  // 面
    bound() = default;
    bound(std::string name_,short type_);
};
inline bound::bound(std::string name_,short type_):name(name_),type(type_) {}

inline std::vector<node_class> NodeList;   // 节点
inline std::vector<cell_class> CellList;   // 网格
inline std::vector<cell_class> GhostList;  // 虚网格
inline std::vector<face_class> FaceList;   // 面
inline std::vector<bound> bounds; // 边界

}

// 节点
namespace cc {
struct node_class {
    int index = 0;                    // 节点编号
    double x = 0.0, y = 0.0, z = 0.0; // 节点坐标
    node_class() = default;
    node_class(int number_, double x_, double y_, double z_);
};
inline node_class::node_class(int number_,double x_,double y_,double z_):index(number_),x(x_),y(y_),z(z_){}
} 

// 面
namespace cc {
struct face_class {
    int index = 0;                // 面编号
    short type = INTER;           // 面类型
    node_class* node[4] = {};     // 面邻接点指针(按右手定则排序)
    vec3 mid = {};                // 面中点坐标
    vec3 nor = {};                // 带面积法向量
    bool outer = false;           // 面法向与nei[0]->nei[1]同向时为true
    double area = 0.0;            // 面面积
    cell_class* nei[2] = {};      // 面邻接网格指针
    vecp phy;                     // 物理量
    vecgrad phgrad;               // 面法向梯度
    struct otphy otphy;           // 引申物理量
    double un = 0.0;              // 面法向速度
    #if SAchosen
        SA::turbulence tur;       // 湍流
    #elif SSTchosen
        SST::turbulence tur;      // 湍流
    #endif
    vecp phynei[2];               // 左右插值物理量
    struct otphy otnei[2];        // 左右引申物理量
    vec5 convect;                 // 无粘通量
    vec5 visflux;                 // 黏性通量

    face_class() = default;
    face_class(int index_, int p1_, int p2_, int p3_, int p4_, int c1_, int c2_, short type_);
};

inline face_class::face_class(int index_, int p1_, int p2_, int p3_, int p4_, int c1_, int c2_, short type_):
    index(index_),type(type_){
        node[0]=&NodeList[p1_-1];node[1]=&NodeList[p2_-1];node[2]=&NodeList[p3_-1];node[3]=&NodeList[p4_-1];
        nei[0] = c1_ > 0 ? &CellList[c1_-1] : nullptr;
        nei[1] = c2_ > 0 ? &CellList[c2_-1] : nullptr;}
}

// 格子
namespace cc {

struct cell_class {
    int index = 0;                             // 编号
    int ecnt = 6;                              // 面邻接边个数
    bool fnorm[6] = {};                        // 邻接面外法向标记
    cell_class* nei[6] = {};                   // 邻接网格指针
    face_class* faces[6] = {};                 // 邻接面指针
    double vol = 0.0;                          // 体积
    vec3 center;                               // 中心坐标
    vec3 proj;                                 // 各方向投影面积和,用于当地时间步长
    double localdt = 0.0;                      // 当地时间步长

    vecp phy;                                   // 物理量
    struct otphy otphy;                         // 引申物理量
    vec5 conser;                                // 守恒量
    vec5 conserformer;                          // 前期守恒量
    vec5 lastconser;                            // 上一时间步守恒量
    vec5 convect;                               // 无粘对流项
    vec5 visflux;                               // 黏性通量
    vecgrad phgrad;                             // 物理量梯度
    #if JSTchosen
        JST::dissipation diss;     // 耗散项
    #elif ROEchosen
        ROE::dissipation diss;     // 耗散项
    #endif
    #if SAchosen
        SA::turbulence tur;       // 湍流
    #elif SSTchosen
        SST::turbulence tur;      // 湍流
    #endif 
    LSCBmatrix LSCB;      // LSCB梯度预处理矩阵

    cell_class() = default;
    cell_class(int index_, int f1_, int f2_, int f3_, int f4_, int f5_, int f6_);
};

inline cell_class* othercell(cell_class* c,face_class* f){
    return f->nei[0] == c ? f->nei[1] : f->nei[0];
}

inline cell_class::cell_class(int index_, int f1_, int f2_, int f3_, int f4_, int f5_, int f6_):index(index_)
    {faces[0]=&FaceList[f1_-1];faces[1]=&FaceList[f2_-1];faces[2]=&FaceList[f3_-1];faces[3]=&FaceList[f4_-1];
    faces[4]=&FaceList[f5_-1];faces[5]=&FaceList[f6_-1];}

} 