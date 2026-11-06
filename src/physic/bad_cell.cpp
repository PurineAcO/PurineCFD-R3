#include "physic.h"

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
        if (nbad < 100) printf("  bad cell #%d (%.6g, %.6g, %.6g)\n", c.index, c.center.x,
                               c.center.y, c.center.z);
        ++nbad;
    }
    if (nbad > 100) printf("  ... 共 %d 个坏单元\n", nbad);
    return nbad;
}