#include "implicit/Implicit.h"
#include "physics/Physics.h"
#include "turbulence/Turbulence.h"

cc::mat5 SA::diffusion_state_jac(const cc::face_class& face){
    const double len = face.nor.norm();
    const double nx = face.nor.x/len, ny = face.nor.y/len, nz = face.nor.z/len;
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const double nu = face.tur->miubl, T = face.phy.T, mu = face.otphy->mu;
    const cc::vecgrad& g = face.phgrad;
    const double div = g.ugrad.x + g.vgrad.y + g.wgrad.z;
    const double sxx = 2.0*g.ugrad.x - (2.0/3.0)*div;
    const double syy = 2.0*g.vgrad.y - (2.0/3.0)*div;
    const double szz = 2.0*g.wgrad.z - (2.0/3.0)*div;
    const double sxy = g.ugrad.y + g.vgrad.x;
    const double sxz = g.ugrad.z + g.wgrad.x;
    const double syz = g.vgrad.z + g.wgrad.y;
    const double fxs = nx*sxx + ny*sxy + nz*sxz;
    const double fys = nx*sxy + ny*syy + nz*syz;
    const double fzs = nx*sxz + ny*syz + nz*szz;
    const double Tn = nx*g.Tgrad.x + ny*g.Tgrad.y + nz*g.Tgrad.z;
    cc::vec5 cs[5];
    prim_sens(face, cs);
    const double chi = _chi(rho, mu, nu);
    const double fv1 = _fv1(chi);
    const double D = fv1 + 3.0*fv1*(1.0 - fv1);
    const double MUT = ph::_MUT(mu, T);
    const cc::vec5 e_rho(1.0, 0.0, 0.0, 0.0, 0.0);
    const cc::vec5 dmueff = MUT*cs[3] + (nu*D)*e_rho;
    const cc::vec5 dmuT = MUT*cs[3];
    const cc::vec5 dmut = dmueff - dmuT;
    const cc::vec5 dk = ph::Cp*((1.0/ph::Pr)*dmuT + (1.0/Prt)*dmut);
    const cc::vec5 r1 = (len*fxs)*dmueff;
    const cc::vec5 r2 = (len*fys)*dmueff;
    const cc::vec5 r3 = (len*fzs)*dmueff;
    const cc::vec5 r4 = (len*fxs*u)*dmueff + (len*fys*v)*dmueff + (len*fzs*w)*dmueff + (len*Tn)*dk;
    cc::mat5 b;
    const cc::vec5* rows[5] = {nullptr, &r1, &r2, &r3, &r4};
    for(int r = 0; r < 5; r++){
        if(rows[r] == nullptr){for(int c = 0; c < 5; c++) b.a[r][c] = 0.0; continue;}
        b.a[r][0] = (*rows[r]).c; b.a[r][1] = (*rows[r]).x; b.a[r][2] = (*rows[r]).y;
        b.a[r][3] = (*rows[r]).z; b.a[r][4] = (*rows[r]).e;
    }
    return b;
}
