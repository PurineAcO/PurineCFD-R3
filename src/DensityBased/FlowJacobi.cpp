#include "turbulence/SA.h"
#include "implict.h"
#include "euler/euler.h"
#include "physic/physic.h"
#include "boundary/FAR.h"
#include "boundary/NoSlipWall.h"
#include "boundary/SYM.h"

void impl::form_euler_jacobi(cc::cell_class& cell){
    cell.conj = form_pfgpw_convect(cell);
}

void impl::form_viscous_jacobi(cc::face_class& face){
    face.imats.Bf = SA::diffusion_state_jac(face);
    if(face.type != cc::INTER){
        face.imats.Ab.clear();
        face.imats.Tf.clear();
        return;
    }
    face.imats.Ab = roe_abs_n(face);
    face.imats.Tf = SA::diffusion_grad_jac(face, SA::normal_dist(&face));
}

void impl::wall_jacobi(cc::face_class* face, cc::cell_class* cell, double s){
    const double len = face->nor.norm();
    const double un = std::abs(face->phy.u*face->nor.x + face->phy.v*face->nor.y +
                               face->phy.w*face->nor.z)/len;
    const double lamc = (un + face->otphy->a)*len;
    cell->diag += 0.5*lamc*cc::eyes;
    cell->diag += SA::diffusion_grad_jac(*face, cell->tur->sad);
    cell->diag += (-s)*face->imats.Bf;
}

void impl::sym_jacobi(cc::face_class* face, cc::cell_class* cell, double s){
    const double len = face->nor.norm();
    const double un = std::abs(face->phy.u*face->nor.x + face->phy.v*face->nor.y +
                               face->phy.w*face->nor.z)/len;
    const double lamc = (un + face->otphy->a)*len;
    cell->diag += 0.5*lamc*cc::eyes;
    cell->diag += SA::diffusion_grad_jac(*face, SA::normal_dist(face));
    cell->diag += (-s)*face->imats.Bf;
}

void impl::far_jacobi(cc::face_class* face, cc::cell_class* cell, double s){
    const double len = face->nor.norm();
    const double un = std::abs(face->phy.u*face->nor.x + face->phy.v*face->nor.y +
                               face->phy.w*face->nor.z)/len;
    const double lamc = (un + face->otphy->a)*len;
    cell->diag += 0.5*lamc*cc::eyes;
    cell->diag += (-s)*face->imats.Bf;
}

double impl::face_lambda(const cc::face_class* f, const cc::cell_class& cell){
    const double len = f->nor.norm();
    const double un = std::abs(f->phy.u*f->nor.x + f->phy.v*f->nor.y + f->phy.w*f->nor.z)/len;
    const double Cn = std::max(4.0/3.0, ph::gamma/ph::Pr);
    return (un + f->otphy->a)*len + Cn*f->tur->mueff/f->phy.rho*len*len/cell.vol;
}

void impl::form_cell_lambda(cc::cell_class& cell){
    cell.lambda = 0.0;
    allface(cell) {cell.lambda += face_lambda(cell.faces[i], cell);}
}

void impl::form_diag(cc::cell_class& cell, double inv_dt, double cfl){
    form_cell_lambda(cell);
    cell.diag.clear();
    cell.diag += (cell.vol*inv_dt + cell.lambda/cfl)*cc::eyes;
    cc::mat5 far;
    allface(cell){
        const double s = 2*cell.fnorm[i] - 1;
        cc::face_class* f = cell.faces[i];
        if(f->type == cc::INTER) {cell.diag += 0.5*s*f->toface_jacobi(cell.conj) + 0.5*f->imats.Ab
                                              + f->imats.Tf + (-0.5*s)*f->imats.Bf;}
        else if(f->type == cc::WALL){wall_jacobi(f, &cell, s);}
        else if(f->type == cc::SYM){sym_jacobi(f, &cell, s);}
        else if(f->type == cc::FAR){far_jacobi(f, &cell, s); far += far_convect_jacobi(f);}
    }
    cell.diag += far;
}