#pragma once
#include <string>
#include <vector>

#include "udf/udf.h"

#define allface(cell) for (int i = 0; i < (cell).fcnt; i++)
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

extern std::vector<node_class> NodeList;   // 节点
extern std::vector<cell_class> CellList;   // 网格
extern std::vector<face_class> FaceList;   // 面
extern std::vector<bound> bounds;          // 边界

inline int cell_num = 0;  // 网格数目
inline int face_num = 0;  // 面数目
inline int node_num = 0;  // 节点数目
inline int type_num = 0;  // 边界类型数目

// 对于边界面找到内部网格
cell_class* boundary_findcell(face_class* face);
// 找到邻接网格
cell_class* othercell(cc::cell_class* c,cc::face_class* f);

}

// 分配指针内存
void allocate_fields();
// 释放指针内存
void release_fields();