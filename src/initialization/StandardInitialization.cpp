#include "initialization/Initialization.h"
#include "data_structure/DataStructure.h"
#include "physics/Physics.h"
#include <cstdio>

void standard_initialization(){
    allocate_fields();
    const cc::vecp farphy = from_boundary_to_phy();
    const double miubl_inf = cc::FAR_config.chi*ph::form_otphy(farphy).mu/farphy.rho;
#pragma omp parallel for schedule(static)
    allcell icell(i).phy = farphy;
#pragma omp parallel for schedule(static)
    allcell *icell(i).otphy = ph::form_otphy(icell(i).phy);
#pragma omp parallel for schedule(static)
    allcell icell(i).tur->miubl = miubl_inf;

    printf("STD Initialization OK!, u is: %f\n",cc::CellList[0].phy.u);

}
