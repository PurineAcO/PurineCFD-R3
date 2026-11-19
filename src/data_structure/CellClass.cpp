#include "data_structure/DataStructure.h"
#include "physics/Physics.h"

void cc::cell_class::prim() {
    phy.rho = conser.c;
    phy.u = conser.x / conser.c;
    phy.v = conser.y / conser.c;
    phy.w = conser.z / conser.c;
    const double e = conser.e / conser.c;
    phy.T = (e - 0.5 * (phy.u * phy.u + phy.v * phy.v + phy.w * phy.w)) / ph::Cv;
    *otphy = ph::form_otphy(phy);
}

void cc::cell_class::form_conservative() {
    conser.c = phy.rho;
    conser.x = phy.rho * phy.u;
    conser.y = phy.rho * phy.v;
    conser.z = phy.rho * phy.w;
    conser.e = phy.rho * ph::get_energy(phy);
    lastconser = conser;
}

cc::cell_class::cell_class(int index_, int f1_, int f2_, int f3_, int f4_,
                           int f5_, int f6_, int fcnt_):
    index(index_), fcnt(fcnt_){
        const int ids[6] = {f1_, f2_, f3_, f4_, f5_, f6_};
        for(int j = 0; j < fcnt_ && j < 6; j++){
            faces[j] = ids[j] > 0 ? &FaceList[ids[j]-1] : nullptr;
        }
}
