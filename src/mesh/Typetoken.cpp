#include "mesh/Mesh.h"
#include "data_structure/DataStructure.h"
#include <cstring>

void typetoken(char *a, char *b){
    if(strcmp(b,"WALL") == 0){cc::bounds.emplace_back(std::string(a),cc::WALL);}
    else if(strcmp(b,"VIL") == 0){cc::bounds.emplace_back(std::string(a),cc::VIL);}
    else if(strcmp(b,"FAR") == 0){cc::bounds.emplace_back(std::string(a),cc::FAR);}
    else if(strcmp(b,"SYM") == 0){cc::bounds.emplace_back(std::string(a),cc::SYM);}
    else if(strcmp(b,"POL") == 0){cc::bounds.emplace_back(std::string(a),cc::POL);}
    else if(strcmp(b,"INTER") == 0){cc::bounds.emplace_back(std::string(a),cc::INTER);}
}
