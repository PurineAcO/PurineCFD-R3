#pragma once
#include "data_structure/face_class.h"
#include "data_structure/cell_class.h"
#include "data_structure/vectors.h"

struct roe_waves{
    cc::vec3 n;         // 界面法向量
    double alpha[5];    // 波强度
    double lambda[5];   // 特征值
    cc::vec5 r[5];      // 右特征向量
};

// 汇总单元各面的无粘和黏性通量
void assemble_flux(cc::cell_class& cell);
// Roe无粘通量
void convect_ROE(cc::face_class& face);
// Roe分解
roe_waves roe_decompose(const cc::face_class& face);
// 形成三向无粘通量Jacobi
cc::mat5tri form_pfgpw_convect(const cc::cell_class& cell);
// Roe通量Jacobian绝对值
cc::mat5 roe_abs_n(const cc::face_class& face);

