#include "LSCB.h"
#include "turbulence/SA.h"
#include "boundary/NoSlipWall.h"

void grad_onface(cc::cell_class& cell){
    allface(cell){
        cc::face_class* face = cell.faces[i];
        if(face->type == cc::INTER){
            const cc::cell_class* lo = face->nei[0];
            const cc::cell_class* hi = face->nei[1];
            const double dlo = cc::dot(face->nor, face->mid - lo->center);
            const double dhi = cc::dot(face->nor, hi->center - face->mid);
            const double w = dlo/(dhi + dlo);
            face->phgrad = w*lo->phgrad + (1.0 - w)*hi->phgrad;
            face->tur->miublgrad = w*lo->tur->miublgrad + (1.0 - w)*hi->tur->miublgrad;
        }else if(face->type == cc::WALL){
            const double y1 = cell.tur->sad;
            cc::vec3 n = face->mid - cell.center;
            n = (1.0/n.norm())*n;
            const cc::vec3 seed = std::abs(n.x) < 0.9 ? cc::vec3(1.0,0.0,0.0) : cc::vec3(0.0,1.0,0.0);
            cc::vec3 t1 = cc::cross(n, seed);
            t1 = (1.0/t1.norm())*t1;
            const cc::vec3 t2 = cc::cross(n, t1);
            auto split = [&](const cc::vec3& g, double wall, double inner){
                const double dn = (wall - inner)/y1;
                return dn*n + cc::dot(g, t1)*t1 + cc::dot(g, t2)*t2;
            };
            face->phgrad.clear();
            face->phgrad.ugrad = split(cell.phgrad.ugrad, face->phy.u, cell.phy.u);
            face->phgrad.vgrad = split(cell.phgrad.vgrad, face->phy.v, cell.phy.v);
            face->phgrad.wgrad = split(cell.phgrad.wgrad, face->phy.w, cell.phy.w);
            face->phgrad.Tgrad = split(cell.phgrad.Tgrad, face->phy.T, cell.phy.T);
            face->tur->miublgrad = split(cell.tur->miublgrad, face->tur->miubl, cell.tur->miubl);
        }else{
            face->phgrad = cell.phgrad;
            face->tur->miublgrad = cell.tur->miublgrad;
        }
    }
}