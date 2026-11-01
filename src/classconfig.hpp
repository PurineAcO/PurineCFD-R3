#include <vector>


namespace cc {

struct cell_class;     // 3D单元
struct face_class;     // 3D面
struct node_class;     // 节点
struct bound;          // 边界
inline std::vector<cell_class> CellList;  // 网格
inline std::vector<face_class> FaceList;  // 面
inline std::vector<node_class> NodeList;  // 节点

inline constexpr short FAR = 0;
inline constexpr short VIL = 1;
inline constexpr short POL = 2;
inline constexpr short INTER = 3;

struct bound{
    const short type;             // 类型
    std::vector<face_class*> fs;  // 面
};
inline std::vector<bound> bounds;  // 边界

} // namespace cc

namespace cc {

struct cell_class{
    int index;              // 编号
    short fcnt;             // 邻接面个数
    node_class* nodes[8];   // 顶点指针
    face_class* faces[6];   // 邻接面指针
    
};


}