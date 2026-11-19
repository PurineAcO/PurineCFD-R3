#include "implicit/Implicit.h"
#include "physics/Physics.h"
#include "turbulence/Turbulence.h"

void SA::source_jacobi_diag_SA(cc::cell_class& cell, double inv_dt, double cfl){    cell.tur->diag = cell.vol*cell.phy.rho*inv_dt;
    cell.tur->diag += cell.phy.rho*cell.lambda/cfl;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0)*f->convect.c;
        cell.tur->diag += std::max(mdot, 0.0);
        cell.tur->diag += Gamma_eff(cell, f)*f->nor.norm()/normal_dist(f);
    }
    const double chidft2 = -2.0*Ct4*cell.tur->chi*cell.tur->chi*cell.tur->ft2;
    const double chidfv2 = -cell.tur->chi*(1.0 - 3.0*cell.tur->chi*cell.tur->fv1*(1.0 - cell.tur->fv1))/
                           std::pow(1.0 + cell.tur->chi*cell.tur->fv1, 2);
    const double cw3_6 = std::pow(Cw3, 6);
    const double g6 = std::pow(cell.tur->g, 6);
    const double h = std::pow((1.0 + cw3_6)/(cw3_6 + g6), 1.0/6.0);
    const double rdfw = cell.tur->r*h*cw3_6/(cw3_6 + g6)*(1.0 + Cw2*(6.0*std::pow(cell.tur->r, 5) - 1.0));
    const double rhoS = cell.phy.rho*cell.tur->Sbl;
    cell.tur->diag -= cell.vol*(rhoS*Cb1*(1.0 - cell.tur->ft2 - chidft2) +
                                rhoS*Cb1*cell.tur->r*((cell.tur->fv2 + chidfv2)*(1.0 - cell.tur->ft2) +
                                                      2.0*cell.tur->ft2 + chidft2) -
                                rhoS*Cw1*kappa*kappa*cell.tur->r*
                                    (2.0*cell.tur->fw + rdfw - cell.tur->r*rdfw*(cell.tur->fv2 + chidfv2)));
    cell.tur->diag += cell.vol*2.0*C5*cell.phy.rho*cell.tur->miubl*cell.tur->S2/
                      (ph::gamma*ph::R*cell.phy.T);
}
