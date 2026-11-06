#pragma once

#include "SA.hpp"
#include "boundary.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "convect.hpp"
#include "dissipation.hpp"
#include "grad.hpp"
#include "interpolate.hpp"
#include <cmath>
#include <cstddef>
#include <algorithm>

#define forrk(rktable) for (int z = 0; z < (int)std::size(rktable); z++)
#define allface_rk(f) for (cc::face_class& f : cc::FaceList)
#define allcell_rk(c) for (cc::cell_class& c : cc::CellList)

namespace RK {
inline constexpr double RK5[5] = {0.25, 1.0 / 6, 0.375, 0.5, 1.0};
inline constexpr double RK3[3] = {0.2075, 0.5915, 1.0};
}

namespace CFLsup {
inline constexpr double Cv = 1.0;
inline const double Cn = (1.333 > cc::gamma / cc::Pr) ? 1.333 : cc::gamma / cc::Pr;
}

// 当地时间步长
void local_timestep(cc::cell_class& cell);
// 稳态形式RK求解1次伪时间
void one_rans(double dt, bool urans);

inline void local_timestep(cc::cell_class& cell){
    double lambda_c = 0.0;
    double lambda_v = 0.0;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const double len = f->nor.norm();
        const double un = f->toface_jacobi(f->phy.u, f->phy.v, f->phy.w)*(2*cell.fnorm[i] - 1)/len;
        lambda_c += (std::abs(un) + f->otphy.a)*len;
        lambda_v += CFLsup::Cv*f->tur.mueff/f->phy.rho*CFLsup::Cn*len*len/cell.vol;
    }
    cell.localdt = fatime::CFL*cell.vol/(lambda_c + lambda_v);
}

inline void one_rans(double dt, bool urans){
    allcell_rk(c) c.copyconver();
    forrk(RK::RK3){
        allcell_rk(c) c.prim();
        noslip_wall_boundary();
        symmetry_wall_boundary();
        far_field_boundary();
        allface_rk(f) interpolate_mid(&f);
        allface_rk(f) muscl(&f);
        allcell_rk(c) least_square_cell_based(c);
        allcell_rk(c) grad_onface(c);
        if(cc::scheme == 'J'){
            allcell_rk(c) jst::shockwave_recognize(c);
            allcell_rk(c) jst::laplace_dissipation(c);
            allface_rk(f) convect_JST(f);
        }else{
            allface_rk(f) convect_ROE(f);
        }
        allcell_rk(c) assemble_flux(c, cc::scheme);
        allface_rk(f) SA::diffusion_SA(f);
        allcell_rk(c) SA::assemble_visflux(c);
        allcell_rk(c) local_timestep(c);
        if(!urans){
            allcell_rk(c) c.conser = c.conserformer -
                RK::RK3[z]/c.vol*c.localdt*(c.convect - c.visflux);
        }else{
            allcell_rk(c) c.conser = c.conserformer -
                RK::RK3[z]/(1.0/c.localdt + 1.0/(2*dt))*
                (1.0/c.vol*(c.convect - c.visflux) + 1.0/(2*dt)*(c.conser - c.lastconser));
        }
    }
    allcell_rk(c) SA::SA_equation_after(c, c.localdt, dt, urans);
}
