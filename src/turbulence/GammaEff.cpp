#include "turbulence/Turbulence.h"

double SA::Gamma_eff(cc::cell_class& c, cc::face_class* f){
    if(f->type == cc::WALL) {return inv_sigma*f->otphy->mu;}
    if(f->type == cc::FAR) {return inv_sigma*(f->otphy->mu + f->phy.rho*f->tur->miubl);}
    if(f->type == cc::INTER){
        return (1.0 + Cb2)*inv_sigma*(f->otphy->mu + f->phy.rho*f->tur->miubl) -
               Cb2*inv_sigma*(c.otphy->mu + c.phy.rho*c.tur->miubl);
    }
    return inv_sigma*(f->otphy->mu + f->phy.rho*f->tur->miubl);
}
