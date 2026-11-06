#pragma once

#include "SA.hpp"
#include "boundary.hpp"
#include "classconfig.hpp"
#include "config.hpp"
#include "convect.hpp"
#include "grad.hpp"
#include "interpolate.hpp"
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace impl {
inline std::vector<cc::vec5> dQ, dQold, rhs;
inline int n_clip = 0;
inline int sweeps = 2;

// 面法向无粘Jacobian(含面积)
cc::mat5 euler_jacobi_n(const cc::face_class& f);
// 面无粘谱半径
double face_radius(const cc::face_class& f);
// 面粘性谱半径
double face_visc_radius(const cc::face_class& f);
// 组装右端残差
void right_residual();
// 分配解向量
void malloc_vector();
// 复制解向量
void Qcopy();
// 形成对角块
void form_diag(cc::cell_class& cell, double inv_dt, double cfl);
// 对称高斯赛德尔单点核
void solvekernel(int c, int lo, int hi);
// 线性系统松弛
int relax();
// 隐式推进一个伪时间步
double flow_step(double dt, double cfl);
// 隐式推进湍流方程
void turbulence_step(double dt, double cfl);
// 导出隐式Jacobian左端算子
bool export_jacobi(const char* path);

template <class KERNEL> void block_sgs(KERNEL&& kernel);

inline cc::mat5 euler_jacobi_n(const cc::face_class& f){
    const double len = f.nor.norm();
    const cc::vec3 n = (1.0/len)*f.nor;
    const double u = f.phy.u, v = f.phy.v, w = f.phy.w;
    const double q2 = u*u + v*v + w*w;
    const double H = cc::Cp*f.phy.T + 0.5*q2;
    const double un = u*n.x + v*n.y + w*n.z;
    const double g1 = cc::gamma - 1.0;
    cc::mat5 A;
    A[0][1] = n.x; A[0][2] = n.y; A[0][3] = n.z;
    A[1][0] = n.x*0.5*g1*q2 - u*un;
    A[1][1] = un + (2.0 - cc::gamma)*u*n.x;
    A[1][2] = u*n.y - g1*v*n.x;
    A[1][3] = u*n.z - g1*w*n.x;
    A[1][4] = g1*n.x;
    A[2][0] = n.y*0.5*g1*q2 - v*un;
    A[2][1] = v*n.x - g1*u*n.y;
    A[2][2] = un + (2.0 - cc::gamma)*v*n.y;
    A[2][3] = v*n.z - g1*w*n.y;
    A[2][4] = g1*n.y;
    A[3][0] = n.z*0.5*g1*q2 - w*un;
    A[3][1] = w*n.x - g1*u*n.z;
    A[3][2] = w*n.y - g1*v*n.z;
    A[3][3] = un + (2.0 - cc::gamma)*w*n.z;
    A[3][4] = g1*n.z;
    A[4][0] = un*(0.5*g1*q2 - H);
    A[4][1] = n.x*H - g1*u*un;
    A[4][2] = n.y*H - g1*v*un;
    A[4][3] = n.z*H - g1*w*un;
    A[4][4] = cc::gamma*un;
    return len*A;
}

inline double face_radius(const cc::face_class& f){
    const double len = f.nor.norm();
    const double un = f.toface_jacobi(f.phy.u, f.phy.v, f.phy.w)/len;
    return (std::abs(un) + f.otphy.a)*len;
}

inline double face_visc_radius(const cc::face_class& f){
    const double len = f.nor.norm();
    const double Cn = std::max(4.0/3.0, cc::gamma/cc::Pr);
    return Cn*f.tur.mueff/f.phy.rho*len*len;
}

inline void right_residual(){
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

inline void malloc_vector(){
    if((int)dQ.size() == cc::cell_num) {return;}
    dQ.resize(cc::cell_num);
    dQold.resize(cc::cell_num);
    rhs.resize(cc::cell_num);
}

inline void Qcopy(){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < (int)dQ.size(); i++) {dQold[i] = dQ[i];}
}

inline void form_diag(cc::cell_class& cell, double inv_dt, double cfl){
    cell.diag.clear();
    double lambda = 0.0;
    allface(cell){
        cc::face_class* f = cell.faces[i];
        const double s = (f->nei[0] == &cell) ? 1.0 : -1.0;
        lambda += face_radius(*f) + face_visc_radius(*f)/cell.vol;
        const double d = 0.5*(face_radius(*f) + face_visc_radius(*f)/cell.vol);
        for(int r = 0; r < 5; r++) {cell.diag[r][r] += d;}
    }
    const double d0 = cell.vol*inv_dt + lambda/cfl;
    for(int r = 0; r < 5; r++) {cell.diag[r][r] += d0;}
}

inline void solvekernel(int c, int lo, int hi){
    const cc::cell_class& cell = cc::CellList[c];
    cc::vec5 acc = rhs[c];
    allface(cell){
        const cc::face_class* f = cell.faces[i];
        if(f->type != cc::INTER) {continue;}
        const cc::cell_class* o = othercell(const_cast<cc::cell_class*>(&cell),
                                            const_cast<cc::face_class*>(f));
        const int j = o->index - 1;
        const cc::vec5 dj = (j >= lo && j < hi) ? dQ[j] : dQold[j];
        const double s = (f->nei[0] == &cell) ? 1.0 : -1.0;
        cc::mat5 off = (0.5*s)*euler_jacobi_n(*f);
        const double d = 0.5*(face_radius(*f) + face_visc_radius(*f)/cell.vol);
        for(int r = 0; r < 5; r++) {off[r][r] -= d;}
        acc -= off*dj;
    }
    cc::mat5 d = cell.diag;
    if(!d.solve(acc)) {acc.clear();}
    dQ[c] = acc;
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

inline int relax(){
    int clipped = 0;
#pragma omp parallel for schedule(static) reduction(+ : clipped)
    for(int i = 0; i < cc::cell_num; i++){
        cc::cell_class& c = cc::CellList[i];
        cc::vec5 q = c.conser + dQ[i];
        const double ke = 0.5*(q.x*q.x + q.y*q.y + q.z*q.z)/q.c;
        const double e = q.e - ke;
        if(!(q.c > 0.0) || !(e > 0.0) || !std::isfinite(e)){
            dQ[i] *= 0.5;
            clipped++;
            continue;
        }
        c.conser = q;
    }
    n_clip += clipped;
    return clipped;
}

inline double flow_step(double dt, double cfl){
    malloc_vector();
    right_residual();
    const double inv_dt = dt > 0.0 ? 1.0/(2.0*dt) : 0.0;
#pragma omp parallel for schedule(static)
    for(int i = 0; i < cc::cell_num; i++){
        const cc::cell_class& c = cc::CellList[i];
        rhs[i] = (-1.0)*(c.convect - c.visflux);
        if(dt > 0.0) {rhs[i] -= (inv_dt*0.5*c.vol)*(c.conser - c.lastconser);}
        dQ[i].clear();
    }
    allcell form_diag(icell(i), inv_dt, cfl);
    for(int sw = 0; sw < sweeps; sw++){Qcopy(); block_sgs(solvekernel);}
    relax();
    double sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : sum)
    for(int i = 0; i < cc::cell_num; i++){
        const cc::vec5& r = rhs[i];
        sum += r.c*r.c + r.x*r.x + r.y*r.y + r.z*r.z + r.e*r.e;
    }
    allcell icell(i).lastconser = icell(i).conser;
    return std::sqrt(sum/cc::cell_num);
}

inline void turbulence_step(double dt, double cfl){
#pragma omp parallel for schedule(static)
    for(int i = 0; i < cc::cell_num; i++){
        cc::cell_class& c = cc::CellList[i];
        double lambda = 0.0;
        allface(c){
            cc::face_class* f = c.faces[i];
            lambda += (face_radius(*f) + face_visc_radius(*f)/c.vol)*SA::inv_sigma;
        }
        c.conserformer = c.conser;
        SA::SA_equation_after(c, cfl*c.vol/lambda, dt, cc::urans);
    }
}

inline bool export_jacobi(const char* path){
    std::vector<long long> row_ptr(cc::cell_num + 1, 0);
    for(int c = 0; c < cc::cell_num; c++){
        row_ptr[c + 1] = row_ptr[c] + 1 + cc::CellList[c].fcnt;
    }
    const long long nnz = row_ptr[cc::cell_num];
    std::vector<long long> col_idx(nnz, 0);
    std::vector<double> vals((size_t)nnz*25, 0.0);
    for(int c = 0; c < cc::cell_num; c++){
        const cc::cell_class& cell = cc::CellList[c];
        long long k = row_ptr[c];
        col_idx[k] = c;
        for(int r = 0; r < 5; r++){
            for(int s = 0; s < 5; s++) {vals[(size_t)k*25 + r*5 + s] = cell.diag[r][s];}
        }
        k++;
        allface(cell){
            const cc::face_class* f = cell.faces[i];
            if(f->type != cc::INTER) {continue;}
            const cc::cell_class* o = othercell(const_cast<cc::cell_class*>(&cell),
                                               const_cast<cc::face_class*>(f));
            col_idx[k] = o->index - 1;
            const double s = (f->nei[0] == &cell) ? 1.0 : -1.0;
            cc::mat5 off = (0.5*s)*euler_jacobi_n(*f);
            const double d = 0.5*(face_radius(*f) + face_visc_radius(*f)/cell.vol);
            for(int r = 0; r < 5; r++) {off[r][r] -= d;}
            for(int r = 0; r < 5; r++){
                for(int q = 0; q < 5; q++) {vals[(size_t)k*25 + r*5 + q] = off[r][q];}
            }
            k++;
        }
    }
    FILE* fp = nullptr;
    fopen_s(&fp, path, "wb");
    if(fp == nullptr) {return false;}
    const long long hdr[2] = {cc::cell_num, nnz};
    fwrite(hdr, sizeof(long long), 2, fp);
    fwrite(row_ptr.data(), sizeof(long long), row_ptr.size(), fp);
    fwrite(col_idx.data(), sizeof(long long), col_idx.size(), fp);
    fwrite(vals.data(), sizeof(double), vals.size(), fp);
    const bool ok = ferror(fp) == 0;
    fclose(fp);
    return ok;
}
}
