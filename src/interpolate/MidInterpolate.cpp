#include "interpolate.h"
#include "physic/physic.h"

void interpolate_mid(cc::face_class* face){
    if(face->type != cc::INTER) {return;}
    face->phy = 0.5*(face->nei[0]->phy + face->nei[1]->phy);

    #if SAchosen
        face->tur->miubl = 0.5*(face->nei[0]->tur->miubl + face->nei[1]->tur->miubl);
    #elif SSTchosen
        /* waiting for SST */
    #endif
    
    *face->otphy = ph::form_otphy(face->phy);
    const cc::vec3 v{face->phy.u, face->phy.v, face->phy.w};
    face->otphy->un = cc::dot(v, face->nor)/face->nor.norm();
}