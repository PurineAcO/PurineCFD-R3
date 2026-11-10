#pragma once

#include "SA.hpp"
#include "boundary.hpp"
#include "convect.hpp"
#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "grad/LSCB.h"
#include "interpolate.hpp"
#include "udf/udf.h"
#include <algorithm>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace impl {
inline std::vector<cc::vec5> rhs, dQ, dQold;
inline std::vector<double> lam, diagv, resid, dnu, dnuold;
inline int n_clip = 0;
inline int sweeps = 2;

// 计算整场残差 R = convect - visflux (单元净流出通量)
void eval_residual();
// 流动方程一次隐式伪时间步, 返回连续方程残差的L2
double flow_step(double dt, double cfl);
// SA方程一次隐式伪时间步, 需紧跟flow_step调用(复用其面通量与梯度)
void sa_step(double dt, double cfl);
// 面的SA扩散系数
double sa_coef(const cc::face_class* f);
// 单元中心到面另一侧沿面法向的距离
double normal_dist(const cc::face_class* f, const cc::cell_class& cell, const cc::cell_class* nb);
// 内部面另一侧的单元, 边界面返回nullptr
cc::cell_class* other_cell(const cc::face_class* f, const cc::cell_class& cell);
// 面谱半径
double face_lambda(const cc::face_class* f, const cc::cell_class& cell);
// 守恒量沿n的无粘通量
cc::vec5 euler_flux(const cc::vec5& Q, double nx, double ny, double nz);
// SA破坏项对nu的导数
double sa_destruction_jac(const cc::cell_class& cell);

inline void eval_residual(){
    allcell icell(i).prim();
    noslip_wall_boundary();
    symmetry_wall_boundary();
    far_field_boundary();
    allfac {interpolate_mid(&iface(i)); muscl(&iface(i));}
    allcell least_square_cell_based(icell(i));
    allcell grad_onface(icell(i));
    allfac convect_ROE(iface(i));
    allcell assemble_flux(icell(i), 'R');
    allfac SA::diffusion_SA(iface(i));
    allcell SA::assemble_visflux(icell(i));
}

inline cc::vec5 euler_flux(const cc::vec5& Q, double nx, double ny, double nz){
    const double u = Q.x/Q.c, v = Q.y/Q.c, w = Q.z/Q.c;
    const double p = (ph::gamma - 1.0)*(Q.e - 0.5*Q.c*(u*u + v*v + w*w));
    const double un = u*nx + v*ny + w*nz;
    return cc::vec5(Q.c*un, Q.x*un + p*nx, Q.y*un + p*ny, Q.z*un + p*nz, (Q.e + p)*un);
}

inline double face_lambda(const cc::face_class* f, const cc::cell_class& cell){
    const double len = f->nor.norm();
    const double un = std::abs(f->phy.u*f->nor.x + f->phy.v*f->nor.y + f->phy.w*f->nor.z)/len;
    const double Cn = std::max(4.0/3.0, ph::gamma/ph::Pr);
    return (un + f->otphy->a)*len + Cn*f->tur->mueff/f->phy.rho*len*len/cell.vol;
}

inline cc::cell_class* other_cell(const cc::face_class* f, const cc::cell_class& cell){
    if(f->type != cc::INTER) {return nullptr;}
    return f->nei[0] == &cell ? f->nei[1] : f->nei[0];
}

inline double sa_coef(const cc::face_class* f){
    return SA::inv_sigma*(f->otphy->mu + f->phy.rho*std::max(f->tur->miubl, 0.0));
}

inline double normal_dist(const cc::face_class* f, const cc::cell_class& cell,
                          const cc::cell_class* nb){
    const cc::vec3 d = nb ? nb->center - cell.center : f->mid - cell.center;
    return std::abs(cc::dot(d, f->nor))/f->nor.norm();
}

template <class KERNEL> inline void block_sgs(KERNEL&& kernel){
#pragma omp parallel
    {
#ifdef _OPENMP
        const int nt = omp_get_num_threads(), tid = omp_get_thread_num();
#else
        const int nt = 1, tid = 0;
#endif
        const int lo = (int)((long long)cc::cell_num*tid/nt);
        const int hi = (int)((long long)cc::cell_num*(tid + 1)/nt);
        for(int c = lo; c < hi; c++) {kernel(c, lo, hi);}
        for(int c = hi - 1; c >= lo; c--) {kernel(c, lo, hi);}
    }
}

inline double flow_step(double dt, double cfl){
    const int N = cc::cell_num;
    if((int)rhs.size() != N){
        rhs.resize(N); dQ.resize(N); dQold.resize(N);
        lam.resize(N); diagv.resize(N); resid.resize(N);
        dnu.resize(N); dnuold.resize(N);
    }
    eval_residual();
    const double inv_dt = dt > 0.0 ? 1.0/dt : 0.0;
    double res2 = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : res2)
    for(int c = 0; c < N; c++){
        cc::cell_class& cell = cc::CellList[c];
        double sum = 0.0;
        allface(cell) {sum += face_lambda(cell.faces[i], cell);}
        lam[c] = sum;
        diagv[c] = cell.vol*inv_dt + sum/cfl + 0.5*sum;
        rhs[c] = cell.convect - cell.visflux;
        if(dt > 0.0) {rhs[c] += (cell.vol*inv_dt)*(cell.conser - cell.lastconser);}
        dQ[c].clear();
        res2 += rhs[c].c*rhs[c].c;
    }
    auto kernel = [](int c, int lo, int hi){
        const cc::cell_class& cell = cc::CellList[c];
        cc::vec5 acc = rhs[c];
        allface(cell){
            const cc::face_class* f = cell.faces[i];
            const cc::cell_class* nb = other_cell(f, cell);
            if(nb == nullptr) {continue;}
            const int j = nb->index - 1;
            const cc::vec5& dj = (j >= lo && j < hi) ? dQ[j] : dQold[j];
            const double s = 2*cell.fnorm[i] - 1;
            const double nx = s*f->nor.x, ny = s*f->nor.y, nz = s*f->nor.z;
            acc += 0.5*(euler_flux(nb->conser + dj, nx, ny, nz) - euler_flux(nb->conser, nx, ny, nz)) -
                   (0.5*face_lambda(f, cell))*dj;
        }
        dQ[c] = acc*(-1.0/diagv[c]);
    };
    for(int sw = 0; sw < sweeps; sw++){
#pragma omp parallel for schedule(static)
        for(int c = 0; c < N; c++) {dQold[c] = dQ[c];}
        block_sgs(kernel);
    }
    int clipped = 0;
#pragma omp parallel for schedule(static) reduction(+ : clipped)
    for(int c = 0; c < N; c++){
        const cc::vec5 Q = cc::CellList[c].conser;
        const double ein = Q.e - 0.5*(Q.x*Q.x + Q.y*Q.y + Q.z*Q.z)/Q.c;
        double w = 1.0;
        for(int it = 0; it < 6; it++){
            const cc::vec5 Qn = Q + w*dQ[c];
            const double ken = 0.5*(Qn.x*Qn.x + Qn.y*Qn.y + Qn.z*Qn.z);
            if(Qn.c > 0.5*Q.c && Qn.e - ken/Qn.c > 0.5*ein) {break;}
            w *= 0.5;
            clipped++;
        }
        cc::CellList[c].conser = Q + w*dQ[c];
    }
    n_clip += clipped;
    return std::sqrt(res2/N);
}

inline double sa_destruction_jac(const cc::cell_class& cell){
    const double nu = std::max(cell.tur->miubl, 0.0);
    const double chi = cell.phy.rho*nu/cell.otphy->mu;
    const double fv1 = chi*chi*chi/(chi*chi*chi + SA::Cv1*SA::Cv1*SA::Cv1);
    const double fv2 = 1.0 - chi/(1.0 + chi*fv1);
    const double d2 = cell.tur->sad*cell.tur->sad;
    const cc::vecgrad& g = cell.phgrad;
    const double cx = g.wgrad.y - g.vgrad.z;
    const double cy = g.ugrad.z - g.wgrad.x;
    const double cz = g.vgrad.x - g.ugrad.y;
    const double vort = std::sqrt(cx*cx + cy*cy + cz*cz);
    const double sn = nu/(d2*SA::kappa*SA::kappa);
    const double St = std::max(vort + fv2*sn, std::max(0.3*vort, 1e-20));
    const double r = std::min(sn/St, SA::rmax);
    const double gg = r + SA::Cw2*(std::pow(r, 6) - r);
    const double c6 = std::pow(SA::Cw3, 6);
    const double fw = gg*std::pow((1.0 + c6)/(std::pow(gg, 6) + c6), 1.0/6.0);
    return 2.0*cell.phy.rho*SA::Cw1*fw*nu/d2;
}

inline void sa_step(double dt, double cfl){
    const int N = cc::cell_num;
    const double inv_dt = dt > 0.0 ? 1.0/dt : 0.0;
#pragma omp parallel for schedule(static)
    for(int c = 0; c < N; c++){
        const cc::cell_class& cell = cc::CellList[c];
        const double nu = cell.tur->miubl;
        double res = 0.0, dg = 0.0;
        allface(cell){
            const cc::face_class* f = cell.faces[i];
            const cc::cell_class* nb = other_cell(f, cell);
            const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0)*f->convect.c;
            const double nu_o = nb ? nb->tur->miubl : f->tur->miubl;
            res += mdot*(mdot > 0.0 ? nu : nu_o);
            const double coef = sa_coef(f);
            double dist = 0.0, dnudn = 0.0;
            if(nb){
                const cc::vec3 d = nb->center - cell.center;
                const double len = d.norm();
                const cc::vec3 e = (1.0/len)*d;
                const cc::vec3 g = f->tur->miublgrad;
                const double corr = (nu_o - nu)/len - cc::dot(g, e);
                dnudn = cc::dot(g + corr*e, (2*cell.fnorm[i] - 1)*f->nor);
                dist = normal_dist(f, cell, nb);
            }else{
                dist = (f->type == cc::WALL) ? cell.tur->sad : normal_dist(f, cell, nullptr);
                dnudn = (nu_o - nu)/dist*f->nor.norm();
            }
            res -= coef*dnudn;
            dg += std::max(mdot, 0.0) + coef*f->nor.norm()/dist;
        }
        res -= cell.vol*SA::source_SA(cell);
        if(dt > 0.0){
            res += cell.vol*inv_dt*(cell.conser.c*nu - cell.lastconser.c*cell.tur->miubl_former);
        }
        dg += cell.conser.c*(cell.vol*inv_dt + lam[c]/cfl) + cell.vol*sa_destruction_jac(cell);
        resid[c] = res;
        diagv[c] = dg;
        dnu[c] = 0.0;
    }
    auto kernel = [](int c, int lo, int hi){
        const cc::cell_class& cell = cc::CellList[c];
        double acc = resid[c];
        allface(cell){
            const cc::face_class* f = cell.faces[i];
            const cc::cell_class* nb = other_cell(f, cell);
            if(nb == nullptr) {continue;}
            const int j = nb->index - 1;
            const double dj = (j >= lo && j < hi) ? dnu[j] : dnuold[j];
            const double mdot = (f->nei[0] == &cell ? 1.0 : -1.0)*f->convect.c;
            acc += (std::min(mdot, 0.0) - sa_coef(f)*f->nor.norm()/normal_dist(f, cell, nb))*dj;
        }
        dnu[c] = -acc/diagv[c];
    };
    for(int sw = 0; sw < sweeps; sw++){
#pragma omp parallel for schedule(static)
        for(int c = 0; c < N; c++) {dnuold[c] = dnu[c];}
        block_sgs(kernel);
    }
#pragma omp parallel for schedule(static)
    for(int c = 0; c < N; c++){
        cc::CellList[c].tur->miubl = std::max(0.0, cc::CellList[c].tur->miubl + dnu[c]);
    }
}

}
