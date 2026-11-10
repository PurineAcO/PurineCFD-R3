#include "LSCB.h"

void cc::LSCBmatrix::form_invA(){
    invA[0][0] = Ayy*Azz - Ayz*Ayz;
    invA[1][1] = Axx*Azz - Axz*Axz;
    invA[2][2] = Axx*Ayy - Axy*Axy;
    invA[0][1] = invA[1][0] = Axz*Ayz - Axy*Azz;
    invA[0][2] = invA[2][0] = Axy*Ayz - Axz*Ayy;
    invA[1][2] = invA[2][1] = Axy*Axz - Ayz*Axx;
    const double det = Axx*invA[0][0] + Axy*invA[0][1] + Axz*invA[0][2];
    for(int r=0;r<3;r++)
        for(int s=0;s<3;s++) invA[r][s] /= det;
}