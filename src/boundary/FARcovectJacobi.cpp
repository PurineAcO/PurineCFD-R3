#include "FAR.h"
#include "physic/physic.h"

cc::mat5 far_convect_jacobi(const cc::face_class* face){
    const double len = face->nor.norm();
    double nx = face->nor.x/len, ny = face->nor.y/len, nz = face->nor.z/len;
    const cc::cell_class* cb = cc::boundary_findcell(const_cast<cc::face_class*>(face));
    if(nx*(face->mid.x - cb->center.x) + ny*(face->mid.y - cb->center.y) +
       nz*(face->mid.z - cb->center.z) < 0.0) {nx = -nx; ny = -ny; nz = -nz;}
    const double rho = face->phy.rho, u = face->phy.u, v = face->phy.v, w = face->phy.w;
    const double p = face->otphy->p;
    const double a = face->otphy->a;
    const double g = ph::gamma, g1 = g - 1.0;
    const double q2 = u*u + v*v + w*w;
    const double H = ph::Cp*face->phy.T + 0.5*q2;
    const double un = u*nx + v*ny + w*nz;
    const bool inflow = un < 0.0;
    const cc::vec3 nv(nx, ny, nz);
    const cc::vec3 seed = std::abs(nx) < 0.9 ? cc::vec3(1.0,0.0,0.0) : cc::vec3(0.0,1.0,0.0);
    cc::vec3 t1 = cc::cross(nv, seed);
    t1 = (1.0/t1.norm())*t1;
    const cc::vec3 t2 = cc::cross(nv, t1);
    const double ut = u*t1.x + v*t1.y + w*t1.z;
    const double wt = u*t2.x + v*t2.y + w*t2.z;
    double A[5][5] = {};
    A[0][0] = un; A[0][1] = rho;
    A[1][0] = u*un;
    A[1][1] = rho*(2.0*nx*un + t1.x*ut + t2.x*wt);
    A[1][2] = rho*t1.x*un; A[1][3] = rho*t2.x*un; A[1][4] = nx;
    A[2][0] = v*un;
    A[2][1] = rho*(2.0*ny*un + t1.y*ut + t2.y*wt);
    A[2][2] = rho*t1.y*un; A[2][3] = rho*t2.y*un; A[2][4] = ny;
    A[3][0] = w*un;
    A[3][1] = rho*(2.0*nz*un + t1.z*ut + t2.z*wt);
    A[3][2] = rho*t1.z*un; A[3][3] = rho*t2.z*un; A[3][4] = nz;
    A[4][0] = 0.5*un*(un*un + ut*ut + wt*wt);
    A[4][1] = rho*(H + un*un);
    A[4][2] = rho*un*ut;
    A[4][3] = rho*un*wt;
    A[4][4] = g/g1*un;
    const double dpQ[5] = {0.5*g1*q2, -g1*u, -g1*v, -g1*w, g1};
    const double duQ[5] = {-u/rho, 1.0/rho, 0.0, 0.0, 0.0};
    const double dvQ[5] = {-v/rho, 0.0, 1.0/rho, 0.0, 0.0};
    const double dwQ[5] = {-w/rho, 0.0, 0.0, 1.0/rho, 0.0};
    double B[5][5] = {};
    for(int k = 0; k < 5; k++){
        const double dvn = nx*duQ[k] + ny*dvQ[k] + nz*dwQ[k];
        const double dvt1 = t1.x*duQ[k] + t1.y*dvQ[k] + t1.z*dwQ[k];
        const double dvt2 = t2.x*duQ[k] + t2.y*dvQ[k] + t2.z*dwQ[k];
        const double da = a/(2.0*p)*dpQ[k] - (k == 0 ? a/(2.0*rho) : 0.0);
        B[0][k] = rho/(2.0*a)*dvn + rho/(g1*a)*da;
        B[1][k] = 0.5*dvn + da/g1;
        B[2][k] = inflow ? 0.0 : dvt1;
        B[3][k] = inflow ? 0.0 : dvt2;
        B[4][k] = g*p/rho*B[0][k];
    }
    if(!inflow){
        const double s_ = p/std::pow(rho, g);
        for(int k = 0; k < 5; k++){
            const double ds = s_/p*dpQ[k] - (k == 0 ? g*s_/rho : 0.0);
            B[0][k] += -rho/(g1*s_)*ds;
            B[4][k] += p/s_*ds;
        }
    }
    cc::mat5 res;
    for(int r = 0; r < 5; r++){
        for(int c = 0; c < 5; c++){
            double s = 0.0;
            for(int m = 0; m < 5; m++) s += A[r][m]*B[m][c];
            res.a[r][c] = s*len;
        }
    }
    return res;
}