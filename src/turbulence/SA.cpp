#include "SA.h"
#include "physic/physic.h"

void SA::diffusion_SA(cc::face_class& face){
    SA::build_SA(face);
    const double mueff = face.phy.rho*face.tur->fv1*face.tur->miubl + face.otphy->mu;
    face.tur->mueff = mueff;
    const cc::vecgrad& g = face.phgrad;
    const double div = g.ugrad.x + g.vgrad.y + g.wgrad.z;
    const double tauxx = mueff*(2.0*g.ugrad.x - (2.0/3.0)*div);
    const double tyy = mueff*(2.0*g.vgrad.y - (2.0/3.0)*div);
    const double tzz = mueff*(2.0*g.wgrad.z - (2.0/3.0)*div);
    const double txy = mueff*(g.ugrad.y + g.vgrad.x);
    const double txz = mueff*(g.ugrad.z + g.wgrad.x);
    const double tyz = mueff*(g.vgrad.z + g.wgrad.y);
    const double mut = face.phy.rho*face.tur->fv1*face.tur->miubl;
    const double lam = face.otphy->mu/ph::Pr + mut/SA::Prt;
    const cc::vec3 q = (-lam*ph::Cp)*g.Tgrad;
    const double u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const cc::vec5 F(0.0, tauxx, txy, txz, u*tauxx + v*txy + w*txz - q.x);
    const cc::vec5 G(0.0, txy, tyy, tyz, u*txy + v*tyy + w*tyz - q.y);
    const cc::vec5 K(0.0, txz, tyz, tzz, u*txz + v*tyz + w*tzz - q.z);
    face.visflux = face.nor.x*F + face.nor.y*G + face.nor.z*K;
}

void SA::assemble_visflux(cc::cell_class& cell){
    cell.visflux.clear();
    allface(cell) {cell.visflux += (2*cell.fnorm[i] - 1)*cell.faces[i]->visflux;}
}

double SA::source_SA(cc::cell_class& cell){
    build_SA(cell);
    const double d2 = cell.tur->sad*cell.tur->sad;
    const double production = Cb1*(1.0 - cell.tur->ft2)*cell.tur->Sbl*cell.phy.rho*cell.tur->miubl;
    const double destruction = cell.phy.rho*(Cw1*cell.tur->fw - Cb1/(kappa*kappa)*cell.tur->ft2)*
                               cell.tur->miubl*cell.tur->miubl*(1.0/d2);
    const double gradient_source = Cb2*inv_sigma*cell.phy.rho*
                                   cc::dot(cell.tur->miublgrad, cell.tur->miublgrad);
    const double compressible = C5*cell.phy.rho*cell.tur->miubl*cell.tur->miubl*cell.tur->S2/
                                (ph::gamma*ph::R*cell.phy.T);
    return production - destruction + gradient_source - compressible;
}