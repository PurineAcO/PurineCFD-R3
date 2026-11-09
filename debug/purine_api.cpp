#include "data_structure/classconfig.h"
#include "data_structure/cell_class.h"
#include "data_structure/face_class.h"
#include "data_structure/node_class.h"
#include "geometry/geometry.h"
#include "initialize/initialization.h"
#include "mesh/readmesh.h"
#include "physic/physic.h"

#ifdef _WIN32
#define PG_API extern "C" __declspec(dllexport)
#else
#define PG_API extern "C" __attribute__((visibility("default")))
#endif

static void reset_mesh() {
    release_fields();
    cc::NodeList.clear();
    cc::FaceList.clear();
    cc::CellList.clear();
    cc::bounds.clear();
    cc::node_num = 0;
    cc::face_num = 0;
    cc::cell_num = 0;
    cc::type_num = 0;
}

PG_API int pg_readmesh(const char* path) {
    reset_mesh();
    return readmesh(path) ? 1 : 0;
}

PG_API void pg_geometry_main() { geometry_main(); }

PG_API int pg_nodes() { return cc::node_num; }

PG_API int pg_faces() { return cc::face_num; }

PG_API int pg_cells() { return cc::cell_num; }

PG_API void pg_copy_nodes(double* xyz) {
    for (int i = 0; i < cc::node_num; i++) {
        xyz[3 * i] = cc::NodeList[i].p.x;
        xyz[3 * i + 1] = cc::NodeList[i].p.y;
        xyz[3 * i + 2] = cc::NodeList[i].p.z;
    }
}

PG_API void pg_copy_face_conn(int* node_id, int* cell_id, int* type, int* pcnt) {
    for (int i = 0; i < cc::face_num; i++) {
        const cc::face_class& f = cc::FaceList[i];
        pcnt[i] = f.pcnt;
        type[i] = f.type;
        for (int k = 0; k < 4; k++) {
            node_id[4 * i + k] = k < f.pcnt ? f.node[k]->index : 0;
        }
        cell_id[2 * i] = f.nei[0] != nullptr ? f.nei[0]->index : 0;
        cell_id[2 * i + 1] = f.nei[1] != nullptr ? f.nei[1]->index : 0;
    }
}

PG_API void pg_copy_cell_conn(int* face_id) {
    for (int i = 0; i < cc::cell_num; i++) {
        const cc::cell_class& c = cc::CellList[i];
        for (int k = 0; k < 6; k++) {
            face_id[6 * i + k] = c.faces[k]->index;
        }
    }
}

PG_API void pg_copy_cell_nei(int* nei_id) {
    for (int i = 0; i < cc::cell_num; i++) {
        const cc::cell_class& c = cc::CellList[i];
        for (int k = 0; k < 6; k++) {
            nei_id[6 * i + k] = c.nei[k] != nullptr ? c.nei[k]->index : 0;
        }
    }
}

PG_API void pg_copy_face_geom(double* area, double* nor, double* mid) {
    for (int i = 0; i < cc::face_num; i++) {
        const cc::face_class& f = cc::FaceList[i];
        area[i] = f.area;
        nor[3 * i] = f.nor.x;
        nor[3 * i + 1] = f.nor.y;
        nor[3 * i + 2] = f.nor.z;
        mid[3 * i] = f.mid.x;
        mid[3 * i + 1] = f.mid.y;
        mid[3 * i + 2] = f.mid.z;
    }
}

PG_API void pg_copy_cell_geom(double* vol, double* cen, int* fnorm) {
    for (int i = 0; i < cc::cell_num; i++) {
        const cc::cell_class& c = cc::CellList[i];
        vol[i] = c.vol;
        cen[3 * i] = c.center.x;
        cen[3 * i + 1] = c.center.y;
        cen[3 * i + 2] = c.center.z;
        for (int k = 0; k < 6; k++) {
            fnorm[6 * i + k] = c.fnorm[k] ? 1 : 0;
        }
    }
}

PG_API double pg_wall_distance(int cell) {
    if (cell < 0 || cell >= cc::cell_num) { return -1.0; }
    return distance_to_wall(cc::CellList[cell]);
}

PG_API void pg_standard_initialization() { standard_initialization(); }

PG_API void pg_from_boundary_to_phy(double* out) {
    const cc::vecp f = from_boundary_to_phy();
    out[0] = f.rho;
    out[1] = f.u;
    out[2] = f.v;
    out[3] = f.w;
    out[4] = f.T;
}

PG_API void pg_far_config(double* out) {
    out[0] = cc::FAR_config.p;
    out[1] = cc::FAR_config.Ma;
    out[2] = cc::FAR_config.T;
    out[3] = cc::FAR_config.AOA;
    out[4] = cc::FAR_config.chi;
}

PG_API void pg_copy_cell_phy(double* phy) {
    for (int i = 0; i < cc::cell_num; i++) {
        const cc::vecp& f = cc::CellList[i].phy;
        phy[5 * i] = f.rho;
        phy[5 * i + 1] = f.u;
        phy[5 * i + 2] = f.v;
        phy[5 * i + 3] = f.w;
        phy[5 * i + 4] = f.T;
    }
}

PG_API void pg_copy_cell_otphy(double* ot) {
    for (int i = 0; i < cc::cell_num; i++) {
        const cc::otphy* o = cc::CellList[i].otphy;
        ot[5 * i] = o->a;
        ot[5 * i + 1] = o->p;
        ot[5 * i + 2] = o->e;
        ot[5 * i + 3] = o->un;
        ot[5 * i + 4] = o->mu;
    }
}

PG_API void pg_copy_cell_miubl(double* miubl) {
    for (int i = 0; i < cc::cell_num; i++) {
        miubl[i] = cc::CellList[i].tur->miubl;
    }
}

PG_API void pg_copy_cell_conser(double* conser) {
    for (int i = 0; i < cc::cell_num; i++) {
        const cc::vec5& q = cc::CellList[i].conser;
        conser[5 * i] = q.c;
        conser[5 * i + 1] = q.x;
        conser[5 * i + 2] = q.y;
        conser[5 * i + 3] = q.z;
        conser[5 * i + 4] = q.e;
    }
}

PG_API void pg_form_conservative() {
    for (int i = 0; i < cc::cell_num; i++) {
        cc::CellList[i].form_conservative();
    }
}

PG_API void pg_prim() {
    for (int i = 0; i < cc::cell_num; i++) {
        cc::CellList[i].prim();
    }
}
