#include "../data_structure/cell_class.h"
#include "../turbulence/SA.h"
#include "../physic/physic.h"
#include "initialization.h"

/*
标准初始化仅仅针对cell,不针对任何face和边界面
*/

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
    // for(cc::face_class& face : cc::FaceList){face.tur->miubl = (face.type == cc::WALL) ? 0.0 : miubl_inf;}
    printf("STD Initialization OK!, u is: %f\n",cc::CellList[0].phy.u);
    
}
