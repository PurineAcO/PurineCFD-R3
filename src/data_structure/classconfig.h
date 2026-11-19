#pragma once
#include <string>
#include <vector>

#include "config/Config.h"
#include "data_structure/Bound.h"

#define allface(cell) for (int i = 0; i < (cell).fcnt; i++)
#define allcell for (int i = 0; i < cc::cell_num; i++)
#define icell(i) cc::CellList[i]
#define allfac for (int i = 0; i < cc::face_num; i++)
#define iface(i) cc::FaceList[i]

namespace cc {

struct cell_class;
struct face_class;
struct node_class;

extern std::vector<node_class> NodeList;
extern std::vector<cell_class> CellList;
extern std::vector<face_class> FaceList;
extern std::vector<bound> bounds;

inline int cell_num = 0;
inline int face_num = 0;
inline int node_num = 0;
inline int type_num = 0;

// 对于边界面找到内部网格
cell_class* boundary_findcell(face_class* face);
// 找到邻接网格
cell_class* othercell(cc::cell_class* c, cc::face_class* f);

}

void malloc_tur(cc::cell_class& c, bool mall);
void malloc_tur(cc::face_class& f, bool mall);
void malloc_otphy(cc::cell_class& c, bool mall);
void malloc_otphy(cc::face_class& f, bool mall);
void malloc_LSCB(cc::cell_class& c, bool mall);
void allocate_fields();
void release_fields();

namespace fatime {
inline double CFL = 1.0;
}
