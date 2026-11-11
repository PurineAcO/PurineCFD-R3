#pragma once

#include "data_structure/vectors.h"
#include "data_structure/classconfig.h"
#include "data_structure/face_class.h"
#include "data_structure/cell_class.h"
#include "data_structure/node_class.h"
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
// 法向量方向校正
double outward_unit_nor(cc::face_class *f, const cc::cell_class *c);
// 几何预处理总流程
void geometry_main();
// 形成一个三维正交基
void tangent_basis(const cc::vec3& n, cc::vec3& t1, cc::vec3& t2);