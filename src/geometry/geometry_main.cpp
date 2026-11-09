#include "geometry.h"

void geometry_main(){
    for(int i = 0; i < cc::cell_num; i++){
        cc::cell_class& c = cc::CellList[i];
        const std::pair<double, cc::vec3> r = cube_volumn_cen(c);
        c.vol = r.first;
        c.center = r.second;
    }
    for(int i = 0; i < cc::face_num; i++){
        cc::face_class& f = cc::FaceList[i];
        f.nor = face_norm(f);
        f.area = face_area(f);
        f.mid = face_mid(f);
    }
    for(int i = 0; i < cc::cell_num; i++){
        cell_face_norm_direction(cc::CellList[i]);
    }
    for(int i = 0; i < cc::cell_num; i++){
        LSCB_pre_process(cc::CellList[i]);
    }
}
