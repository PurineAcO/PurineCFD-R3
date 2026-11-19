#include "implicit/Implicit.h"

bool ifnotneedrelax(const cc::vec5& newed, const cc::vec5& old){
    const double oldein = old.e - 0.5*(old.x*old.x + old.y*old.y + old.z*old.z)/old.c;
    const double newein = newed.e - 0.5*(newed.x*newed.x + newed.y*newed.y + newed.z*newed.z)/newed.c;
    return (newed.c > 0.5*old.c && newein > 0.5*oldein);
}
