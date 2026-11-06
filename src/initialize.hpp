#pragma once
#include "classconfig.hpp"
#include "config.hpp"
#include "physic.hpp"
#include <cmath>

// 标准初始化
void standard_initialization();
// 从流场初始化
void fromflowfield_initialization();
// 从边界条件形成物理量
cc::vecp from_boundary_to_phy();

inline cc::vecp from_boundary_to_phy(){
    const double a = std::sqrt(cc::gamma*cc::R*cc::FAR_config.T);
    const double v = cc::FAR_config.Ma*a;
    return cc::vecp(cc::FAR_config.p/(cc::R*cc::FAR_config.T),
                    v*std::cos(cc::FAR_config.AOA),v*std::sin(cc::FAR_config.AOA),0.0,
                    cc::FAR_config.T);
}

inline void standard_initialization(){
    const cc::vecp farphy = from_boundary_to_phy();
    const double miubl_inf = cc::FAR_config.chi*form_otphy(farphy).mu/farphy.rho;
#pragma omp parallel for schedule(static)
    allcell icell(i).phy = farphy;
#pragma omp parallel for schedule(static)
    allcell icell(i).otphy = form_otphy(icell(i).phy);
#pragma omp parallel for schedule(static)
    allcell icell(i).tur.miubl = miubl_inf;
    for(cc::face_class& face : cc::FaceList){
        face.tur.miubl = (face.type == cc::WALL) ? 0.0 : miubl_inf;
    }
    printf("STD Initialization OK!, u is: %f\n",cc::CellList[0].phy.u);
}

