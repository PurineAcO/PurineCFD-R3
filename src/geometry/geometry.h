#pragma once

#include "../data_structure/vectors.h"
#include "../data_structure/classconfig.h"
#include "../data_structure/face_class.h"
#include "../data_structure/cell_class.h"
#include "../data_structure/node_class.h"
#include <utility>
#include <vector>

// 四面体体积
double tetrahedron_volumn(cc::vec3 a,cc::vec3 b,cc::vec3 c,cc::vec3 d);
// 四面体质心
cc::vec3 tetrahedron_cen(cc::vec3 a,cc::vec3 b,cc::vec3 c,cc::vec3 d);
// 六面体体积和质心
std::pair<double, cc::vec3> cube_volumn_cen(cc::cell_class& c);
// 六面体顶点排序
std::vector<cc::node_class*> sort_cube_node(cc::cell_class& c);
// 面带面积法向
cc::vec3 face_norm(cc::face_class& f);
// 面积
double face_area(cc::face_class &f);
// 面中心
cc::vec3 face_mid(cc::face_class& f);
// 建立面法向关系
void cell_face_norm_direction(cc::cell_class& c);
// 壁面距离
double distance_to_wall(const cc::cell_class& c);
// LSCB预备
void LSCB_pre_process(cc::cell_class& c);




// inline void LSCB_pre_process(cc::cell_class &c){
//     c.LSCB.clear();
//     allface(c){
//         cc::cell_class* nei = cc::othercell(&c, c.faces[i]);
//         if(nei == nullptr) {continue;}
//         const cc::vec3 d = nei->center - c.center;
//         c.LSCB.wi[i] = 1.0/cc::dot(d,d);
//         c.LSCB.Axx += c.LSCB.wi[i] * d.x * d.x;
//         c.LSCB.Axy += c.LSCB.wi[i] * d.x * d.y;
//         c.LSCB.Ayy += c.LSCB.wi[i] * d.y * d.y;
//         c.LSCB.Ayz += c.LSCB.wi[i] * d.y * d.z;
//         c.LSCB.Azz += c.LSCB.wi[i] * d.z * d.z;
//         c.LSCB.Axz += c.LSCB.wi[i] * d.x * d.z;
//     }
//     c.LSCB.form_invA();
// }
