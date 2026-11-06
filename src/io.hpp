#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <cstdio>
#include <string>

// 输出流场
bool dump_field(const char* path);
// 输出壁面压力系数
bool write_wall_cp(const char* path, double pref, const cc::vecp& refp);

inline bool dump_field(const char* path){
    FILE* fp = nullptr;
    fopen_s(&fp, path, "w");
    if(fp == nullptr) {return false;}
    fprintf(fp, "x,y,z,rho,u,v,w,p,T,miubl\n");
    for(const cc::cell_class& c : cc::CellList){
        fprintf(fp, "%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
                c.center.x, c.center.y, c.center.z, c.phy.rho, c.phy.u, c.phy.v, c.phy.w,
                c.otphy.p, c.phy.T, c.tur.miubl);
    }
    const bool ok = ferror(fp) == 0;
    fclose(fp);
    return ok;
}

inline bool write_wall_cp(const char* path, double pref, const cc::vecp& refp){
    FILE* fp = nullptr;
    fopen_s(&fp, path, "w");
    if(fp == nullptr) {return false;}
    const double q = dynamic_pressure(refp);
    fprintf(fp, "x,y,z,theta,cp\n");
    for(const cc::face_class* w : cc::WallFaces){
        const double cp = (w->otphy.p - pref)/q;
        fprintf(fp, "%.9e,%.9e,%.9e,%.6f,%.9e\n", w->mid.x, w->mid.y, w->mid.z,
                std::atan2(w->mid.y, w->mid.x), cp);
    }
    const bool ok = ferror(fp) == 0;
    fclose(fp);
    return ok;
}
