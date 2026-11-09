#include "physic.h"
#include <cstdio>

int bad_cells() {
    int nbad = 0;
    for (const cc::cell_class& c : cc::CellList) {
        const double ke = 0.5 * (c.conser.x * c.conser.x + c.conser.y * c.conser.y +
                                 c.conser.z * c.conser.z) / c.conser.c;
        const double T = (c.conser.e - ke) / (c.conser.c * ph::Cv);
        const bool bad = !std::isfinite(c.conser.c) || c.conser.c <= 0.0 ||
                         !std::isfinite(c.conser.e) || !std::isfinite(c.tur->miubl) ||
                         !(T > 0.0) || !std::isfinite(T);
        if (!bad) { continue; }
        printf("  bad cell #%d (%.6g, %.6g, %.6g)\n", c.index, c.center.x, c.center.y, c.center.z);
        if(c.conser.c <= 0.0) printf("density < 0");
        else if(T <= 0.0) printf("temperature < 0");
        else printf("conservative vector has infinite number");
        ++nbad;
    }
    if (nbad > 100) printf("  ... total: %d bad cells\n", nbad);
    return nbad;
}