#include "face_class.h"

/*
multiply the vector or physics by jacobi matrix
*/

double cc::face_class::toface_jacobi(double F, double G, double K) const {
    return nor.x*F + nor.y*G + nor.z*K;
}

cc::vec5 cc::face_class::toface_jacobi(const cc::vec5& F, const cc::vec5& G, const cc::vec5& K) const {
    return nor.x*F + nor.y*G + nor.z*K;
}
