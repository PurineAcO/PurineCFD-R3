#include "data_structure/classconfig.h"
#include "geometry/geometry.h"
#include "mesh/readmesh.h"
#include <cstdio>

int main(int argc, char** argv){
    const char* meshpath = argc > 1 ? argv[1] : "cases/cyc/mesh.txt";
    if(!readmesh(meshpath)) {return 1;}
    geometry_main();
    std::printf("mesh: %d nodes, %d faces, %d cells\n", cc::node_num, cc::face_num, cc::cell_num);
    return 0;
}
