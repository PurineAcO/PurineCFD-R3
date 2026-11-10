#include "data_structure/vectors.h"
#include "LSCB.h"

cc::vec3 cc::LSCBmatrix::solve(cc::vec3 dphi){
    return cc::vec3(invA[0][0]*dphi.x + invA[0][1]*dphi.y + invA[0][2]*dphi.z,
                    invA[1][0]*dphi.x + invA[1][1]*dphi.y + invA[1][2]*dphi.z,
                    invA[2][0]*dphi.x + invA[2][1]*dphi.y + invA[2][2]*dphi.z);
}