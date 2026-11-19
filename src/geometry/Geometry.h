#pragma once

#include "data_structure/DataStructure.h"

#include <utility>
#include <vector>

// 四面体体积
double tetrahedron_volumn(cc::vec3 a,cc::vec3 b,cc::vec3 c,cc::vec3 d);
// 四面体质心
cc::vec3 tetrahedron_cen(cc::vec3 a,cc::vec3 b,cc::vec3 c,cc::vec3 d);
// 六面体体积和质心
std::pair<double, cc::vec3> cube_volumn_cen(cc::cell_class& c);
#if chose2D
// 三角形面积（二维多边形面积用）
double tri_area_2d(const cc::vec3& a, const cc::vec3& b, const cc::vec3& c);
// 三角形质心（二维多边形形心用）
cc::vec3 tri_cen_2d(const cc::vec3& a, const cc::vec3& b, const cc::vec3& c);
// 多边形面积和质心
std::pair<double, cc::vec3> cube_volumn_cen_2D(cc::cell_class& c);
#endif
// 六面体顶点排序，为使各点满足计算多面体体积所需的对应关系
std::vector<cc::node_class*> sort_cube_node(cc::cell_class& c);
// 带上面积的面法向
// Warning:
// 1. 调用前务必保证四边形（如果提供的是四边形）四点共面，若发生轻微翘曲请先执行修复函数
// 2. 返回的面法向是带面积的单位向量
cc::vec3 face_norm(cc::face_class& f);
#if !chose2D
// 三角形面积（面面积积分用）
double triangle_area(cc::vec3 a,cc::vec3 b,cc::vec3 c);
#endif
// 面积
double face_area(cc::face_class &f);
#if !chose2D
// 三角形面积（面中心加权用）
double ta(cc::vec3 a,cc::vec3 b,cc::vec3 c);
#endif
// 面中心
cc::vec3 face_mid(cc::face_class& f);
// 建立面法向关系
// fnorm 为真表示 cell.faces[i] 的法向指向单元外
// 使用 cell.fnorm[i] 时应当用 (2*fnorm[i]-1) 标记外法向
void cell_face_norm_direction(cc::cell_class& c);
// 壁面距离，启用湍流模型时需要；采用全局搜索求最近壁面距离
double distance_to_wall(const cc::cell_class& c);
// 法向量方向校正：无法从 face 快速转到 cell，所以需要本函数
double outward_unit_nor(cc::face_class *f, const cc::cell_class *c);
// 几何预处理总流程
void geometry_main();
// 形成一个三维正交基
void tangent_basis(const cc::vec3& n, cc::vec3& t1, cc::vec3& t2);
