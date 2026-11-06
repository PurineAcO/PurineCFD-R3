#pragma once

#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
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

// 形成单元湍流中间变量
void build_SA(cc::cell_class& c);
// 形成面湍流中间变量
void build_SA(cc::face_class& f);
// 形成面上粘性通量
void diffusion_SA(cc::face_class& face);
// 组装单元粘性通量
void assemble_visflux(cc::cell_class& cell);
// 形成SA源项
double source_SA(cc::cell_class& cell);
// 推进湍流输运方程
void SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans);
// 面上有效扩散系数
double Gamma_eff(cc::cell_class& c, cc::face_class* f);

inline void build_SA(cc::cell_class& c){
    const cc::vecgrad& g = c.phgrad;
    c.otphy.mu = sutherland::dynamic_viscosity(c.phy.T);
    c.tur.chi = c.phy.rho*(c.tur.miubl > 0.0 ? c.tur.miubl : 0.0)/c.otphy.mu;
    const double chi3 = c.tur.chi*c.tur.chi*c.tur.chi;
    c.tur.fv1 = chi3/(chi3 + Cv1*Cv1*Cv1);
    c.tur.ft2 = Ct3*std::exp(-Ct4*c.tur.chi*c.tur.chi);
    c.tur.fv2 = 1.0 - c.tur.chi/(1.0 + c.tur.chi*c.tur.fv1);
    const double wx = g.vgrad.z - g.wgrad.y;
    const double wy = g.wgrad.x - g.ugrad.z;
    const double wz = g.ugrad.y - g.vgrad.x;
    c.tur.Omega = std::sqrt(wx*wx + wy*wy + wz*wz);
    c.tur.Sbl = std::max(c.tur.Omega + c.tur.fv2*c.tur.miubl/(c.tur.sad*c.tur.sad)/(kappa*kappa),
                         std::max(0.3*c.tur.Omega, 1e-20));
    c.tur.r = std::min(c.tur.miubl/(c.tur.sad*c.tur.sad)/(kappa*kappa)/c.tur.Sbl, rmax);
    const double r2 = c.tur.r*c.tur.r, r3 = r2*c.tur.r;
    c.tur.g = c.tur.r + Cw2*(r3*r3 - c.tur.r);
    const double cw36 = Cw3*Cw3*Cw3*Cw3*Cw3*Cw3;
    c.tur.fw = c.tur.g*std::pow((1.0 + cw36)/(cw36 + std::pow(c.tur.g, 6)), 1.0/6.0);
    c.tur.S2 = 2.0*(g.ugrad.x*g.ugrad.x + g.vgrad.y*g.vgrad.y + g.wgrad.z*g.wgrad.z) +
               (g.ugrad.y + g.vgrad.x)*(g.ugrad.y + g.vgrad.x) +
               (g.ugrad.z + g.wgrad.x)*(g.ugrad.z + g.wgrad.x) +
               (g.vgrad.z + g.wgrad.y)*(g.vgrad.z + g.wgrad.y);
}

inline void build_SA(cc::face_class& f){
    f.otphy.mu = sutherland::dynamic_viscosity(f.phy.T);
    f.tur.chi = f.phy.rho*(f.tur.miubl > 0.0 ? f.tur.miubl : 0.0)/f.otphy.mu;
    const double chi3 = f.tur.chi*f.tur.chi*f.tur.chi;
    f.tur.fv1 = chi3/(chi3 + Cv1*Cv1*Cv1);
}

inline void diffusion_SA(cc::face_class& face){
    build_SA(face);
    const cc::vecgrad& g = face.phgrad;
    const double mueff = face.phy.rho*face.tur.fv1*face.tur.miubl + face.otphy.mu;
    face.tur.mueff = mueff;
    const double div = g.ugrad.x + g.vgrad.y + g.wgrad.z;
    const double txx = mueff*(2.0*g.ugrad.x - (2.0/3.0)*div);
    const double tyy = mueff*(2.0*g.vgrad.y - (2.0/3.0)*div);
    const double tzz = mueff*(2.0*g.wgrad.z - (2.0/3.0)*div);
    const double txy = mueff*(g.ugrad.y + g.vgrad.x);
    const double txz = mueff*(g.ugrad.z + g.wgrad.x);
    const double tyz = mueff*(g.vgrad.z + g.wgrad.y);
    const double lam = face.otphy.mu/cc::Pr + face.phy.rho*face.tur.fv1*face.tur.miubl/Prt;
    const cc::vec3 q = (-lam*cc::Cp)*g.Tgrad;
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

inline double source_SA(cc::cell_class& cell){
    build_SA(cell);
    const double production = Cb1*(1.0 - cell.tur.ft2)*cell.tur.Sbl*cell.phy.rho*cell.tur.miubl;
    const double destruction = cell.phy.rho*(Cw1*cell.tur.fw - Cb1/(kappa*kappa)*cell.tur.ft2)*
                               cell.tur.miubl*cell.tur.miubl/(cell.tur.sad*cell.tur.sad);
    const double compressible = C5*cell.phy.rho*cell.tur.miubl*cell.tur.miubl*cell.tur.S2/
                                (cc::gamma*cc::R*cell.phy.T);
    return production - destruction - compressible;
}

inline void SA_equation_after(cc::cell_class& cell, double dtau, double dt, bool urans){
    double convect = 0.0;
    double diffusion = 0.0;
    allface(cell){
        cc::face_class* face = cell.faces[i];
        const double mu = sutherland::dynamic_viscosity(face->phy.T);
        const double outer = 2*cell.fnorm[i] - 1;
        const double coef = inv_sigma*(mu + face->phy.rho*face->tur.miubl);
        const double flux = face->phy.rho*face->tur.miubl*
                            (face->phy.u*face->nor.x + face->phy.v*face->nor.y + face->phy.w*face->nor.z);
        convect += outer*flux;
        diffusion += outer*coef*(face->tur.miublgrad.x*face->nor.x +
                                 face->tur.miublgrad.y*face->nor.y +
                                 face->tur.miublgrad.z*face->nor.z);
    }
    const double rhs = (diffusion - convect)/cell.vol + source_SA(cell);
    double new_miubl = 0.0;
    if(!urans){
        new_miubl = (cell.conserformer.c*cell.tur.miubl + relax*dtau*rhs)/cell.phy.rho;
    }else{
        new_miubl = (cell.conserformer.c*cell.tur.miubl +
                     (rhs - (cell.tur.miubl - cell.tur.miubl_former)/(2.0*dt))/
                         (1.0/dtau + 1.0/(2.0*dt)))/cell.phy.rho;
    }
    cell.tur.miubl = std::max(0.0, new_miubl);
}

inline double Gamma_eff(cc::cell_class& c, cc::face_class* f){
    if(f->type == cc::WALL) {return inv_sigma*f->otphy.mu;}
    if(f->type == cc::FAR) {return inv_sigma*(f->otphy.mu + f->phy.rho*f->tur.miubl);}
    if(f->type == cc::INTER){
        return (1.0 + Cb2)*inv_sigma*(f->otphy.mu + f->phy.rho*f->tur.miubl) -
               Cb2*inv_sigma*(c.otphy.mu + c.phy.rho*c.tur.miubl);
    }
    return 0.0;
}

}
