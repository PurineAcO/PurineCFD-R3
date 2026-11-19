#include "implicit/Implicit.h"
#include "euler/Euler.h"

void impl::form_euler_jacobi(cc::cell_class& cell){
    cell.conj = form_pfgpw_convect(cell);
}
