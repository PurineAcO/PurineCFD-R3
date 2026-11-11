#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "implict.h"
#include "turbulence/SA.h"
#include "boundary/FAR.h"
#include "boundary/SYM.h"
#include "boundary/NoSlipWall.h"
#include "interpolate/interpolate.h"
#include "euler/euler.h"
#include "grad/LSCB.h"

double SA::right_residual(cc::cell_class& cell){
    double res = 0.0;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const cc::cell_class* nb = cc::othercell(&cell, f);
        const double nu_o = nb ? nb->tur->miubl : f->tur->miubl;
        const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0)*f->convect.c;
        res += mdot*(mdot > 0.0 ? cell.tur->miubl : nu_o);
        res += Gamma_eff(cell, f)*f->nor.norm()/normal_dist(f)*(cell.tur->miubl - nu_o);
    }
    res -= cell.vol*source_SA(cell);
    return res;
}

void impl::right_residual(){
#pragma omp parallel for schedule(static)
    allcell icell(i).prim();
    noslip_wall_boundary();
    symmetry_wall_boundary();
    far_field_boundary();
#pragma omp parallel for schedule(static)
    allfac {interpolate_mid(&iface(i)); muscl(&iface(i));}
#pragma omp parallel for schedule(static)
    allcell least_square_cell_based(icell(i));
#pragma omp parallel for schedule(static)
    allcell grad_onface(icell(i));
#pragma omp parallel for schedule(static)
    allfac convect_ROE(iface(i));
#pragma omp parallel for schedule(static)
    allcell assemble_flux(icell(i));
#pragma omp parallel for schedule(static)
    allfac SA::diffusion_SA(iface(i));
#pragma omp parallel for schedule(static)
    allcell SA::assemble_visflux(icell(i));
}