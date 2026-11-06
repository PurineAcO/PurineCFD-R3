#pragma once

#include "classconfig.hpp"
#include <cmath>

// GGCB梯度
void green_gauss_cell_based(cc::cell_class& cell);
// LSCB梯度
void least_square_cell_based(cc::cell_class& cell);

inline void green_gauss_cell_based(cc::cell_class& cell){
    cell.phgrad.clear();
    cell.tur.miublgrad = cc::vec3{};
    allface(cell){
        const double s = (2*cell.fnorm[i] - 1)/cell.vol;
        cell.phgrad = cell.phgrad + (cell.faces[i]->phy*s)*cell.faces[i]->nor;
        cell.tur.miublgrad += (cell.faces[i]->tur.miubl*s)*cell.faces[i]->nor;
    }
}

inline void least_square_cell_based(cc::cell_class& cell){
    cc::vec3 brho{}, bu{}, bv{}, bw{}, bT{}, bmiubl{};
    allface(cell){
        cc::cell_class* nei = cell.nei[i];
        if(nei == nullptr) {continue;}
        const cc::vec3 d = (nei->center - cell.center)*cell.LSCB.wi[i];
        brho += (nei->phy.rho - cell.phy.rho)*d;
        bu += (nei->phy.u - cell.phy.u)*d;
        bv += (nei->phy.v - cell.phy.v)*d;
        bw += (nei->phy.w - cell.phy.w)*d;
        bT += (nei->phy.T - cell.phy.T)*d;
        bmiubl += (nei->tur.miubl - cell.tur.miubl)*d;
    }
    cell.phgrad.rhograd = cell.LSCB.solve(brho);
    cell.phgrad.ugrad = cell.LSCB.solve(bu);
    cell.phgrad.vgrad = cell.LSCB.solve(bv);
    cell.phgrad.wgrad = cell.LSCB.solve(bw);
    cell.phgrad.Tgrad = cell.LSCB.solve(bT);
    cell.tur.miublgrad = cell.LSCB.solve(bmiubl);
}
