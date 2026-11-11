#pragma once
#include "data_structure/vectors.h"
#include "data_structure/face_class.h"
#include "data_structure/cell_class.h"
#include "data_structure/node_class.h"


namespace cc {

struct LSCBmatrix {
    double Axx = 0.0;    // 对角xx
    double Axy = 0.0;    // 上三角xy
    double Axz = 0.0;    // 上三角xz
    double Ayy = 0.0;    // 对角yy
    double Ayz = 0.0;    // 上三角yz
    double Azz = 0.0;    // 对角zz
    double wi[6] = {};   // 范数权
    double invA[3][3] = {};   // 几何逆

    LSCBmatrix() = default;
    // 清空
    void clear();
    // 形成几何逆
    void form_invA();
    // 解梯度方程
    cc::vec3 solve(cc::vec3 dphi);
};

inline void LSCBmatrix::clear(){Axx = Ayy = Azz = Axy = Axz = Ayz = 0.0; }

}

// 网格的LSCB梯度
void least_square_cell_based(cc::cell_class& cell);
// 湍流变量的LSCB
void LSCB_turbulence(cc::cell_class& cell);
// 面上退化梯度
void grad_onface(cc::cell_class& cell);


// // GGCB梯度
// void green_gauss_cell_based(cc::cell_class& cell);

// inline void green_gauss_cell_based(cc::cell_class& cell){
//     cell.phgrad.clear();
//     cell.tur.miublgrad = cc::vec3{};
//     allface(cell){
//         const double s = (2*cell.fnorm[i] - 1)/cell.vol;
//         cell.phgrad = cell.phgrad + (cell.faces[i]->phy*s)*cell.faces[i]->nor;
//         cell.tur.miublgrad += (cell.faces[i]->tur.miubl*s)*cell.faces[i]->nor;
//     }
// }


