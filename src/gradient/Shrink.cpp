// 原 src/grad/BarthJespersen.cpp（static shrink 提升为外部链接）
#include "gradient/Gradient.h"

double shrink(double q, double lo, double hi, double dq){
    double psi = 1.0;
    if(dq > 0.0) {psi = (hi - q)/dq;}
    else if(dq < 0.0) {psi = (lo - q)/dq;}
    if(psi > 1.0) {psi = 1.0;}
    if(psi < 0.0) {psi = 0.0;}
    return psi;
}
