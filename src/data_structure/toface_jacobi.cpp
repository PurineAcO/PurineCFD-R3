#include "node_class.h"
#include "cell_class.h"
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

cc::mat5 cc::face_class::toface_jacobi(const mat5tri& m) const {
    return nor.x*m.x + nor.y*m.y + nor.z*m.z;
}