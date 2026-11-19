#include <cmath>

#include "io/Io.h"

bool write_wall_cp_mid(FILE* fp, int step, double pref, const cc::vecp& refp){
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
    int nsel = 0;
    for(const cc::bound& b : cc::bounds){
        if(b.type != cc::WALL) {continue;}
        for(const cc::face_class* w : b.fs){
            if(std::abs(w->mid.z - zmid) > best + 1e-9) {continue;}
            const double cp = (w->otphy->p - pref)/q;
            fprintf(fp, "%d,%.6f,%.9e\n", step, std::atan2(w->mid.y, w->mid.x), cp);
            nsel++;
        }
    }
    static bool reported = false;
    if(!reported){
        reported = true;
        std::printf("  wall_cp_mid: zlo %.6f zhi %.6f zmid %.6f best %.9f selected %d\n",
                    zlo, zhi, zmid, best, nsel);
    }
    fflush(fp);
    return ferror(fp) == 0;
}
