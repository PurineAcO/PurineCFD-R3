#include "mesh/Mesh.h"
#include "data_structure/DataStructure.h"
#include <cstdio>
#include <cstring>

bool readmesh_2D(const char* meshpath){
    printf("the path of mesh: %s\n", meshpath);
    FILE* meshfp = nullptr;
    if(fopen_s(&meshfp, meshpath, "r") != 0 || meshfp == nullptr){printf("not found mesh!\n"); return false;}
    char buf[256];
    if(fscanf_s(meshfp, "%d%d%d%d", &cc::node_num, &cc::face_num, &cc::cell_num, &cc::type_num) != 4){
        fclose(meshfp);
        return false;
    }
    printf("%d %d %d %d\n", cc::node_num, cc::face_num, cc::cell_num, cc::type_num);
    cc::NodeList.resize(cc::node_num);
    cc::FaceList.resize(cc::face_num);
    cc::CellList.resize(cc::cell_num);
    cc::bounds.reserve(cc::type_num);
    for(int i = 0; i < cc::type_num; i++){
        char name_before[30] = {}, name_after[30] = {}, dyh[2] = {};
        if(fscanf_s(meshfp, "%29s %1s %29s", name_before, (unsigned)sizeof(name_before),
                    dyh, (unsigned)sizeof(dyh), name_after, (unsigned)sizeof(name_after)) != 3){
            fclose(meshfp);
            return false;
        }
        printf("%s = %s\n", name_before, name_after);
        typetoken(name_before, name_after);
    }
    while(fgets(buf, sizeof(buf), meshfp) != nullptr){
        if(strncmp(buf, "(node)", 6) == 0) {break;}
    }
    for(int i = 0; i < cc::node_num; i++){
        int index = 0; double x = 0.0, y = 0.0;
        if(fscanf_s(meshfp, "%d %lf %lf", &index, &x, &y) != 3){fclose(meshfp); return false;}
        cc::NodeList[i] = cc::node_class(index, x, y, 0.0);
    }
    while(fgets(buf, sizeof(buf), meshfp) != nullptr){
        if(strncmp(buf, "(edge)", 6) == 0) {break;}
    }
    int facecnt = 0, typecnt = 0;
    while(facecnt < cc::face_num){
        int index = 0, p1 = 0, p2 = 0, c1 = 0, c2 = 0;
        const int got = fscanf_s(meshfp, "%d%d%d%d%d", &index, &p1, &p2, &c1, &c2);
        if(got == 0){
            if(fgets(buf, sizeof(buf), meshfp) == nullptr){fclose(meshfp); return false;}
            if(strncmp(buf, "(end)", 5) == 0){typecnt++; continue;}
            if(typecnt >= (int)cc::bounds.size()){fclose(meshfp); return false;}
            if(strncmp(buf, cc::bounds[typecnt].name.c_str(), cc::bounds[typecnt].name.size()) == 0){continue;}
            fclose(meshfp);
            return false;
        }
        if(got != 5){fclose(meshfp); return false;}
        if(index < 1 || index > cc::face_num){fclose(meshfp); return false;}
        cc::FaceList[index-1] = cc::face_class(index, p1, p2, 0, 0, c1, c2,
                                               cc::bounds[typecnt].type, 2);
        if(cc::bounds[typecnt].type != cc::INTER){
            cc::bounds[typecnt].fs.push_back(&cc::FaceList[index-1]);
        }
        facecnt++;
    }
    while(fgets(buf, sizeof(buf), meshfp) != nullptr){
        if(strncmp(buf, "(cell)", 6) == 0) {break;}
    }
    for(int i = 0; i < cc::cell_num; i++){
        int index = 0, f1 = 0, f2 = 0, f3 = 0, f4 = 0;
        if(fscanf_s(meshfp, "%d%d%d%d%d", &index, &f1, &f2, &f3, &f4) != 5){
            fclose(meshfp);
            return false;
        }
        if(index < 1 || index > cc::cell_num){fclose(meshfp); return false;}
        cc::CellList[index-1] = cc::cell_class(index, f1, f2, f3, f4, 0, 0, 4);
    }
    fclose(meshfp);
    printf("2D mesh: %d nodes, %d faces, %d cells\n", cc::node_num, cc::face_num, cc::cell_num);
    return linkmesh();
}
