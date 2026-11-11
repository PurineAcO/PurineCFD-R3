#include "interpolate.h"
#include "grad/limit.h"
#include "physic/physic.h"

void muscl(cc::face_class* face){
    if(face->type == cc::INTER){
        face->phynei[0] = BarthJespersen(face->nei[0], face->mid - face->nei[0]->center);
        face->phynei[1] = BarthJespersen(face->nei[1], face->mid - face->nei[1]->center);
    }else{
        face->phynei[0] = face->phy;
        face->phynei[1] = face->phy;
    }
    *face->otnei[0] = ph::form_otphy(face->phynei[0]);
    *face->otnei[1] = ph::form_otphy(face->phynei[1]);
}