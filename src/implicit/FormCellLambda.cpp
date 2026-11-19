#include "implicit/Implicit.h"

void impl::form_cell_lambda(cc::cell_class& cell){
    cell.lambda = 0.0;
    allface(cell) {cell.lambda += face_lambda(cell.faces[i], cell);}
}
