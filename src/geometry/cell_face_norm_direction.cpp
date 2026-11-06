#include "geometry.h"

/* 
fnorm is true if the norm vector of cell.faces[i] is towards outside the cell 
when use cell.fnorm[i], you should use (2*fnorm[i]-1) to mark outside 
*/

void cell_face_norm_direction(cc::cell_class &c){
    allface(c) c.fnorm[i] = cc::dot(c.faces[i]->mid-c.center, c.faces[i]->nor) > 0 ? true : false;
}
