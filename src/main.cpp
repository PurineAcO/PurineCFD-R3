// #include <iostream>
#include "geometry.hpp"
#include "readmesh.hpp"

int main(){
    if(!readmesh("cases/cyc/mesh.txt")) return 1;
    geometry_main();
    return 0;
}
