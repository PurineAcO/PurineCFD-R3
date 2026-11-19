#include "viscous/Viscous.h"

void SA::assemble_visflux(cc::cell_class& cell){
    cell.visflux.clear();
    allface(cell) {cell.visflux += (2*cell.fnorm[i] - 1)*cell.faces[i]->visflux;}
}
