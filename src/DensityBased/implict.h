#include "data_structure/cell_class.h"

namespace impl {

inline int n_clip = 0;    // 松弛单元数

// 守恒量变化量
inline std::vector<cc::vec5> dQ;
// 冻结守恒量变化量
inline std::vector<cc::vec5> dQold;
// 右侧残差
inline std::vector<cc::vec5> rhs;

// 形成右侧残差
void right_residual();
// 分配求解器向量
void malloc_vector();
// 清理矩阵
void clear_mat();
// 形成无粘通量Jacobi
void form_euler_jacobi(cc::cell_class& cell);
// 形成面上粘性通量Jacobi
void form_viscous_jacobi(cc::face_class& face);
// 壁面特别处理
void wall_jacobi(cc::face_class* face, cc::cell_class* cell, double s);
// 远场特别处理
void far_jacobi(cc::face_class* face, cc::cell_class* cell, double s);
// 对称面特别处理
void sym_jacobi(cc::face_class* face, cc::cell_class* cell, double s);
// 伪时间谱半径
double face_lambda(const cc::face_class* f, const cc::cell_class& cell);
// 组装单元谱半径
void form_cell_lambda(cc::cell_class& cell);
// 求出对角块
void form_diag(cc::cell_class& cell, double inv_dt, double cfl);
// 复制解向量到缓存
void Qcopy();
// 求解核心
void solvekernel(int c, int lo, int hi);

// 松弛
int relax();
// 进行多次隐式伪时间迭代
double flow_step(double dt, double cfl);
}


namespace SA {
// miubl变化量
inline std::vector<double> dnu;
// 上一步miubl变化量
inline std::vector<double> dnuold;
// miubl残差
inline std::vector<double> nurhs;

// 分配求解器向量
void malloc_vector();
// 单元SA残差
double right_residual(cc::cell_class& cell);
// SA对称高斯赛德尔核
void solvekernel(int c, int lo, int hi);
// 复制解向量到缓存
void nucopy();
// 松弛
void relaxer();
// 进行多次湍流模型伪时间迭代
void turbulence_step(double dt, double cfl);
} // namespace SA
