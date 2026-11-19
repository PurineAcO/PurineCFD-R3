#include "residual/Residual.h"

#include "boundary/Boundary.h"
#include "euler/Euler.h"
#include "gradient/Gradient.h"
#include "interpolation/Interpolation.h"
#include "viscous/Viscous.h"

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
