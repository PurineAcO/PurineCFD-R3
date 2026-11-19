#include "io/Io.h"

bool dump_field(const char* path){
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
