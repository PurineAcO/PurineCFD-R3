#include "geometry/Geometry.h"

double distance_to_wall(const cc::cell_class& c){
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
