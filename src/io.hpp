#pragma once

#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "physic/physic.h"
#include <cmath>
#include <cstdio>

// 输出流场
bool dump_field(const char* path);
// 输出壁面压力系数
bool write_wall_cp(const char* path, double pref, const cc::vecp& refp);
// 输出中截面壁面压力系数(逐步追加)
bool write_wall_cp_mid(FILE* fp, int step, double pref, const cc::vecp& refp);

inline bool dump_field(const char* path){
    FILE* fp = nullptr;
    fopen_s(&fp, path, "w");
    if(fp == nullptr) {return false;}
    fprintf(fp, "x,y,z,rho,u,v,w,p,T,miubl\n");
    for(const cc::cell_class& c : cc::CellList){
        fprintf(fp, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                c.center.x, c.center.y, c.center.z, c.phy.rho, c.phy.u, c.phy.v, c.phy.w,
                c.otphy->p, c.phy.T, c.tur->miubl);
    }
    const bool ok = ferror(fp) == 0;
    fclose(fp);
    return ok;
}

inline bool write_wall_cp(const char* path, double pref, const cc::vecp& refp){
    FILE* fp = nullptr;
    fopen_s(&fp, path, "w");
    if(fp == nullptr) {return false;}
    const double q = ph::dynamic_pressure(refp);
    fprintf(fp, "x,y,z,theta,cp\n");
    for(const cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) {continue;}
        for(const cc::face_class* w : b.fs){
            const double cp = (w->otphy->p - pref)/q;
            fprintf(fp, "%.9e,%.9e,%.9e,%.6f,%.9e\n", w->mid.x, w->mid.y, w->mid.z,
                    std::atan2(w->mid.y, w->mid.x), cp);
        }
    }
    const bool ok = ferror(fp) == 0;
    fclose(fp);
    return ok;
}

inline bool write_wall_cp_mid(FILE* fp, int step, double pref, const cc::vecp& refp){
    if(fp == nullptr) {return false;}
    const double q = ph::dynamic_pressure(refp);
    double zlo = 1e300, zhi = -1e300;
    for(const cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) {continue;}
        for(const cc::face_class* w : b.fs){
            if(w->mid.z < zlo) {zlo = w->mid.z;}
            if(w->mid.z > zhi) {zhi = w->mid.z;}
        }
    }
    const double zmid = 0.5*(zlo + zhi);
    double best = 1e300;
    for(const cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) {continue;}
        for(const cc::face_class* w : b.fs){
            const double d = std::abs(w->mid.z - zmid);
            if(d < best) {best = d;}
        }
    }
    for(const cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) {continue;}
        for(const cc::face_class* w : b.fs){
            if(std::abs(w->mid.z - zmid) > best + 1e-9) {continue;}
            const double cp = (w->otphy->p - pref)/q;
            fprintf(fp, "%d,%.6f,%.9e\n", step, std::atan2(w->mid.y, w->mid.x), cp);
        }
    }
    fflush(fp);
    return ferror(fp) == 0;
}
