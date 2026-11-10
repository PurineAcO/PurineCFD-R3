#pragma once

#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "physic/physic.h"
#include "udf/udf.h"
#include <algorithm>
#include <cmath>

namespace SA {
inline constexpr double Cb1 = 0.1355;
inline constexpr double Cb2 = 0.622;
inline constexpr double Cw1 = 3.2391;
inline constexpr double Cw2 = 0.3;
inline constexpr double Cw3 = 2.0;
inline constexpr double Cv1 = 7.1;
inline constexpr double Ct3 = 1.2;
inline constexpr double Ct4 = 0.5;
inline constexpr double kappa = 0.41;
inline constexpr double inv_sigma = 1.5;
inline constexpr double rmax = 10.0;
inline constexpr double relax = 0.8;
inline constexpr double Prt = 0.9;
inline constexpr double C5 = 3.5;

// 形成面上粘性通量
void diffusion_SA(cc::face_class& face);
// 组装单元粘性通量
void assemble_visflux(cc::cell_class& cell);
// 形成SA源项
double source_SA(const cc::cell_class& cell);
// 推进湍流输运方程
void SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans);
// 湍流粘度比
double _chi(double rho, double mu, double miubl);
// 粘度阻尼函数
double _fv1(double chi);
// 转捩控制函数
double _ft2(double chi);
// 涡量辅助函数
double _fv2(double chi);
// 无量纲距离II
double _g(double r);
// 破坏项控制函数
double _fw(double g);

inline double _chi(double rho, double mu, double miubl){
    return rho*(miubl > 0.0 ? miubl : 0.0)/mu;
}

inline double _fv1(double chi){
    return chi*chi*chi/(chi*chi*chi + Cv1*Cv1*Cv1);
}

inline double _ft2(double chi){
    return Ct3*std::exp(-Ct4*chi*chi);
}

inline double _fv2(double chi){
    return 1.0 - chi/(1.0 + chi*_fv1(chi));
}

inline double _g(double r){
    return r + Cw2*(r*r*r*r*r*r - r);
}

inline double _fw(double g){
    constexpr double cw3_squared = Cw3*Cw3;
    constexpr double cw3_sixth = cw3_squared*cw3_squared*cw3_squared;
    const double g_squared = g*g;
    const double g_sixth = g_squared*g_squared*g_squared;
    return g*std::pow((1.0 + cw3_sixth)/(cw3_sixth + g_sixth), 1.0/6.0);
}

inline void diffusion_SA(cc::face_class& face){
    const double mu = face.otphy->mu;
    const double chi = _chi(face.phy.rho, mu, face.tur->miubl);
    const double mut = face.phy.rho*_fv1(chi)*face.tur->miubl;
    const double mueff = mut + mu;
    face.tur->mueff = mueff;
    const cc::vecgrad& g = face.phgrad;
    const double div = g.ugrad.x + g.vgrad.y + g.wgrad.z;
    const double txx = mueff*(2.0*g.ugrad.x - (2.0/3.0)*div);
    const double tyy = mueff*(2.0*g.vgrad.y - (2.0/3.0)*div);
    const double tzz = mueff*(2.0*g.wgrad.z - (2.0/3.0)*div);
    const double txy = mueff*(g.ugrad.y + g.vgrad.x);
    const double txz = mueff*(g.ugrad.z + g.wgrad.x);
    const double tyz = mueff*(g.vgrad.z + g.wgrad.y);
    const double lam = mu/ph::Pr + mut/Prt;
    const cc::vec3 q = (-lam*ph::Cp)*g.Tgrad;
    const double u = face.phy.u, v = face.phy.v, w = face.phy.w;
    const cc::vec5 F(0.0, txx, txy, txz, u*txx + v*txy + w*txz - q.x);
    const cc::vec5 G(0.0, txy, tyy, tyz, u*txy + v*tyy + w*tyz - q.y);
    const cc::vec5 K(0.0, txz, tyz, tzz, u*txz + v*tyz + w*tzz - q.z);
    face.visflux = face.nor.x*F + face.nor.y*G + face.nor.z*K;
}

inline void assemble_visflux(cc::cell_class& cell){
    cell.visflux.clear();
    allface(cell) {cell.visflux += (2*cell.fnorm[i] - 1)*cell.faces[i]->visflux;}
}

inline double source_SA(const cc::cell_class& cell){
    const double mu = cell.otphy->mu;
    const double chi = _chi(cell.phy.rho, mu, cell.tur->miubl);
    const double ft2 = _ft2(chi);
    const double fv2 = _fv2(chi);
    const cc::vecgrad& g = cell.phgrad;
    const double cx = g.wgrad.y - g.vgrad.z;
    const double cy = g.ugrad.z - g.wgrad.x;
    const double cz = g.vgrad.x - g.ugrad.y;
    const double vorticity = std::sqrt(cx*cx + cy*cy + cz*cz);
    const double d2 = cell.tur->sad*cell.tur->sad;
    const double scaled_nu = cell.tur->miubl*(1.0/d2)/(kappa*kappa);
    const double modified_vorticity =
        std::max(vorticity + fv2*scaled_nu, std::max(0.3*vorticity, 1e-20));
    const double production = Cb1*(1.0 - ft2)*modified_vorticity*cell.phy.rho*cell.tur->miubl;
    const double r = std::min(scaled_nu/modified_vorticity, rmax);
    const double destruction = cell.phy.rho*(Cw1*_fw(_g(r)) - Cb1/kappa/kappa*ft2)*
                               cell.tur->miubl*cell.tur->miubl*(1.0/d2);
    const double gradient_source = Cb2*inv_sigma*cell.phy.rho*
                                   cc::dot(cell.tur->miublgrad, cell.tur->miublgrad);
    const double S2 = 2.0*(g.ugrad.x*g.ugrad.x + g.vgrad.y*g.vgrad.y + g.wgrad.z*g.wgrad.z) +
                      (g.ugrad.y + g.vgrad.x)*(g.ugrad.y + g.vgrad.x) +
                      (g.ugrad.z + g.wgrad.x)*(g.ugrad.z + g.wgrad.x) +
                      (g.vgrad.z + g.wgrad.y)*(g.vgrad.z + g.wgrad.y);
    const double compressible = C5*cell.phy.rho*cell.tur->miubl*cell.tur->miubl*S2/
                                (ph::gamma*ph::R*cell.phy.T);
    return production - destruction + gradient_source - compressible;
}

inline void SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans){
    double convect = 0.0;
    double diffusion = 0.0;
    double divergence = 0.0;
    allface(cell){
        cc::face_class* face = cell.faces[i];
        const double mu = ph::dynamic_viscosity(face->phy.T);
        const double outer = 2*cell.fnorm[i] - 1;
        const double coef = inv_sigma*(mu + face->phy.rho*face->tur->miubl);
        const double vn = face->phy.u*face->nor.x + face->phy.v*face->nor.y +
                          face->phy.w*face->nor.z;
        const double outward_vn = outer*vn;
        double upwind = cell.tur->miubl;
        if(outward_vn < 0.0){
            const cc::cell_class* nb = cell.nei[i];
            upwind = nb ? nb->tur->miubl : face->tur->miubl;
        }
        convect += outward_vn*face->phy.rho*upwind;
        diffusion += outer*coef*(face->tur->miublgrad.x*face->nor.x +
                                 face->tur->miublgrad.y*face->nor.y +
                                 face->tur->miublgrad.z*face->nor.z);
        divergence += outward_vn;
    }
    const double rhs = (diffusion - convect)/cell.vol + source_SA(cell) -
                       cell.phy.rho*cell.tur->miubl*divergence/cell.vol;
    double new_miubl = 0.0;
    if(!urans){
        new_miubl = (cell.conserformer.c*cell.tur->miubl + relax*dtau*rhs)/cell.phy.rho;
    }else{
        new_miubl = (cell.conserformer.c*cell.tur->miubl +
                     (rhs - (cell.tur->miubl - cell.tur->miubl_former)/(2.0*dt))/
                         (1.0/dtau + 1.0/(2.0*dt)))/cell.phy.rho;
    }
    cell.tur->miubl = std::max(0.0, new_miubl);
}

}
