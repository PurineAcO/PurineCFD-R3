#include "LSCB.h"

void least_square_cell_based(cc::cell_class &cell){
    cc::vec3 brho{}, bu{}, bv{}, bw{}, bT{};
    allface(cell){
        cc::cell_class* nei = cell.nei[i];
        if(nei == nullptr) {continue;}
        const cc::vec3 d = (nei->center - cell.center)*cell.LSCB->wi[i];
        brho += (nei->phy.rho - cell.phy.rho)*d;
        bu += (nei->phy.u - cell.phy.u)*d;
        bv += (nei->phy.v - cell.phy.v)*d;
        bw += (nei->phy.w - cell.phy.w)*d;
        bT += (nei->phy.T - cell.phy.T)*d;
    }
    cell.phgrad.rhograd = cell.LSCB->solve(brho);
    cell.phgrad.ugrad = cell.LSCB->solve(bu);
    cell.phgrad.vgrad = cell.LSCB->solve(bv);
    cell.phgrad.wgrad = cell.LSCB->solve(bw);
    cell.phgrad.Tgrad = cell.LSCB->solve(bT);
    LSCB_turbulence(cell);
}