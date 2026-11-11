#include "SA.h"
#include "boundary/NoSlipWall.h"
#include "boundary/FAR.h"
#include "physic/physic.h"

double SA::Gamma_eff(cc::cell_class& c, cc::face_class* f){
    if(f->type == cc::WALL) {return inv_sigma*f->otphy->mu;}
    if(f->type == cc::FAR) {return inv_sigma*(f->otphy->mu + f->phy.rho*f->tur->miubl);}
    if(f->type == cc::INTER){
        return (1.0 + Cb2)*inv_sigma*(f->otphy->mu + f->phy.rho*f->tur->miubl) -
               Cb2*inv_sigma*(c.otphy->mu + c.phy.rho*c.tur->miubl);
    }
    return inv_sigma*(f->otphy->mu + f->phy.rho*f->tur->miubl);
}

double SA::normal_dist(const cc::face_class* f){
    const double len = f->nor.norm();
    if(f->type == cc::INTER && f->nei[0] != nullptr && f->nei[1] != nullptr){
        return std::abs(cc::dot(f->nei[0]->center - f->nei[1]->center, f->nor))/len;
    }
    const cc::cell_class* c = cc::boundary_findcell(const_cast<cc::face_class*>(f));
    if(f->type == cc::WALL) {return c->tur->sad;}
    return std::abs(cc::dot(f->mid - c->center, f->nor))/len;
}

static void prim_sens(const cc::face_class& face, cc::vec5 cs[5]){
    const double rho = face.phy.rho, u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const double T = face.phy.T, nu = face.tur->miubl, Cv = ph::Cv;
    const double q2 = u*u + v*v + w*w;
    cs[0] = cc::vec5(-u/rho, 1.0/rho, 0.0, 0.0, 0.0);
    cs[1] = cc::vec5(-v/rho, 0.0, 1.0/rho, 0.0, 0.0);
    cs[2] = cc::vec5(-w/rho, 0.0, 0.0, 1.0/rho, 0.0);
    cs[3] = cc::vec5(q2/(2.0*Cv*rho) - T/rho, -u/(Cv*rho), -v/(Cv*rho), -w/(Cv*rho), 1.0/(Cv*rho));
    cs[4] = cc::vec5(-nu/rho, 0.0, 0.0, 0.0, 1.0/rho);
}

cc::mat5 SA::diffusion_grad_jac(const cc::face_class& face, double dn){
    const double len = face.nor.norm();
    const double nx = face.nor.x/len, ny = face.nor.y/len, nz = face.nor.z/len;
    const double u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const double mu = face.otphy->mu;
    const double mut = std::max(face.tur->mueff - mu, 0.0);
    const double m = face.tur->mueff;
    const double k = ph::Cp*(mu/ph::Pr + mut/Prt);
    cc::vec5 cs[5];
    prim_sens(face, cs);
    cc::mat5 a;
    for(int i = 0; i < 5; i++){
        const double du = cs[0][i], dv = cs[1][i], dw = cs[2][i];
        const double divi = nx*du + ny*dv + nz*dw;
        const double txx = (2.0*nx*du - (2.0/3.0)*divi)/dn;
        const double tyy = (2.0*ny*dv - (2.0/3.0)*divi)/dn;
        const double tzz = (2.0*nz*dw - (2.0/3.0)*divi)/dn;
        const double txy = (ny*du + nx*dv)/dn;
        const double txz = (nz*du + nx*dw)/dn;
        const double tyz = (nz*dv + ny*dw)/dn;
        const double fx = m*(nx*txx + ny*txy + nz*txz);
        const double fy = m*(nx*txy + ny*tyy + nz*tyz);
        const double fz = m*(nx*txz + ny*tyz + nz*tzz);
        a.a[0][i] = 0.0;
        a.a[1][i] = len*fx;
        a.a[2][i] = len*fy;
        a.a[3][i] = len*fz;
        a.a[4][i] = len*(u*fx + v*fy + w*fz + k*cs[3][i]/dn);
    }
    return a;
}

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