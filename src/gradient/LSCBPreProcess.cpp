// 原 src/grad/LSCB_prepare.cpp
#include "gradient/Gradient.h"

void LSCB_pre_process(cc::cell_class &c){
    c.LSCB->clear();
    allface(c){
        cc::cell_class* nei = cc::othercell(&c, c.faces[i]);
        if(nei == nullptr) {continue;}
        const cc::vec3 d = nei->center - c.center;
        c.LSCB->wi[i] = 1.0/cc::dot(d,d);
        c.LSCB->Axx += c.LSCB->wi[i] * d.x * d.x;
        c.LSCB->Axy += c.LSCB->wi[i] * d.x * d.y;
        c.LSCB->Ayy += c.LSCB->wi[i] * d.y * d.y;
        c.LSCB->Ayz += c.LSCB->wi[i] * d.y * d.z;
        c.LSCB->Azz += c.LSCB->wi[i] * d.z * d.z;
        c.LSCB->Axz += c.LSCB->wi[i] * d.x * d.z;
    }
#if chose2D
    c.LSCB->Azz = 1.0;
    c.LSCB->Ayz = 0.0;
    c.LSCB->Axz = 0.0;
#endif
    c.LSCB->form_invA();
}
