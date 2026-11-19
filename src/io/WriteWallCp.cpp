#include <cmath>

#include "io/Io.h"

bool write_wall_cp(const char* path, double pref, const cc::vecp& refp){
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
