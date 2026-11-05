#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#define triangle_area(a,b,c) (0.5*std::abs(cc::cross((a)-(b),(a)-(c)).norm()))


// 四面体体积
double tetrahedron_volumn(cc::vec3 a,cc::vec3 b,cc::vec3 c,cc::vec3 d);
// 四面体质心
cc::vec3 tetrahedron_cen(cc::vec3 a,cc::vec3 b,cc::vec3 c,cc::vec3 d);
// 六面体体积和质心
std::pair<double, cc::vec3> cube_volumn_cen(cc::cell_class& c);
// 六面体顶点排序
std::vector<cc::node_class*> sort_cube_node(cc::cell_class& c);
// 面单位法向
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
// 几何分析
void geometry_main();

inline double tetrahedron_volumn(cc::vec3 a, cc::vec3 b, cc::vec3 c, cc::vec3 d){
    return std::abs(1.0/6*dot(b-a,cc::cross(c-a, d-a)));
}

inline cc::vec3 tetrahedron_cen(cc::vec3 a, cc::vec3 b, cc::vec3 c, cc::vec3 d){
    return 0.25 * (a+b+c+d);
}

inline std::pair<double, cc::vec3> cube_volumn_cen(cc::cell_class& c){
    const std::vector<cc::node_class*> nd = sort_cube_node(c);
    static constexpr short ndi[6] = {1, 2, 3, 7, 4, 5};
    double vol = 0.0; cc::vec3 cen = {};
    for(int i=0;i<6;i++) {
        double dvol =  tetrahedron_volumn(nd[0]->p, nd[ndi[i]]->p, nd[ndi[(i+1)%6]]->p, nd[6]->p);
        cen += tetrahedron_cen(nd[0]->p, nd[ndi[i]]->p, nd[ndi[(i+1)%6]]->p, nd[6]->p)*dvol;
        vol += dvol;
    }
    return std::pair<double,cc::vec3>(vol,1.0/vol*cen);
}

inline std::vector<cc::node_class*> sort_cube_node(cc::cell_class &c){
    const cc::face_class* bottom = c.faces[0];
    const cc::face_class* top = c.faces[1];
    int shift = 0, flip = 0;
    double best = -1.0;
    for(int f = 0; f < 2; f++){
        for(int s = 0; s < 4; s++){
            double cost = 0.0;
            for(int i = 0; i < 4; i++){
                const int j = f ? (s - i + 8) % 4 : (s + i) % 4;
                const cc::vec3 d = top->node[j]->p - bottom->node[i]->p;
                cost += dot(d, d);
            }
            if(best < 0.0 || cost < best){best = cost; shift = s; flip = f;}
        }
    }
    std::vector<cc::node_class*> node_sorted(8, nullptr);
    for(int i = 0; i < 4; i++){
        node_sorted[i] = bottom->node[i];
        node_sorted[4 + i] = top->node[flip ? (shift - i + 8) % 4 : (shift + i) % 4];
    }
    return node_sorted;
}

inline cc::vec3 face_norm(cc::face_class &f){
    cc::vec3 nor_withs = cc::cross(f.node[0]->p-f.node[1]->p, f.node[0]->p-f.node[2]->p);
    return 1.0/nor_withs.norm()*nor_withs;
}

inline double face_area(cc::face_class &f){
    if(f.pcnt == 3) return triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p);
    return triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p) + 
           triangle_area(f.node[0]->p,f.node[3]->p,f.node[2]->p);
}

inline cc::vec3 face_mid(cc::face_class &f){
    if(f.pcnt == 3) return 1.0/3*(f.node[0]->p+f.node[1]->p+f.node[2]->p);
    else{return 1.0/f.area*
             (1.0/3*(f.node[0]->p+f.node[1]->p+f.node[2]->p)*triangle_area(f.node[0]->p,f.node[1]->p,f.node[2]->p) + 
               1.0/3*(f.node[0]->p+f.node[3]->p+f.node[2]->p)*triangle_area(f.node[0]->p,f.node[3]->p,f.node[2]->p));
    }
}

inline void cell_face_norm_direction(cc::cell_class &c){
    /* fnorm is true if the norm vector of cell.faces[i] is towards outside the cell */
    /* when use cell.fnorm[i], you should use (2*fnorm[i]-1) to mark outside */
    allface(c) c.fnorm[i] = cc::dot(c.faces[i]->mid-c.center, c.faces[i]->nor) > 0 ? true : false;
}

inline double distance_to_wall(const cc::cell_class& c){
    double sad = 1e300;
    for(const cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) continue;
        for(const cc::face_class* f : b.fs){
            const double d = (c.center-f->mid).norm();
            if(d < sad) sad = d;
        }
    }
    return sad;
}

inline void LSCB_pre_process(cc::cell_class &c){
    c.LSCB.clear();
    allface(c){
        cc::cell_class* nei = cc::othercell(&c, c.faces[i]);
        if(nei == nullptr) {continue;}
        c.LSCB.wi[i] = 1.0/std::sqrt((nei->center - c.center).norm());
        c.LSCB.Axx += c.LSCB.wi[i] * (nei->center.x - c.center.x) * (nei->center.x - c.center.x);
        c.LSCB.Axy += c.LSCB.wi[i] * (nei->center.x - c.center.x) * (nei->center.y - c.center.y);
        c.LSCB.Ayy += c.LSCB.wi[i] * (nei->center.y - c.center.y) * (nei->center.y - c.center.y);
        c.LSCB.Ayz += c.LSCB.wi[i] * (nei->center.y - c.center.y) * (nei->center.z - c.center.z);
        c.LSCB.Azz += c.LSCB.wi[i] * (nei->center.z - c.center.z) * (nei->center.z - c.center.z);
        c.LSCB.Axz += c.LSCB.wi[i] * (nei->center.x - c.center.x) * (nei->center.z - c.center.z);
    }
    c.LSCB.form_invA();
}
                                                                                                                                                                         
inline void geometry_main(){
    allcell {
        std::pair<double,cc::vec3> res = cube_volumn_cen(icell(i));
        icell(i).center = res.second;
        icell(i).vol = res.first;
    }
    allfac {
        iface(i).nor = face_norm(iface(i)); 
        iface(i).area = face_area(iface(i));
        iface(i).mid = face_mid(iface(i));
    }
    allcell cell_face_norm_direction(icell(i));
    allcell icell(i).tur.sad = distance_to_wall(icell(i));
    allcell LSCB_pre_process(icell(i));
}