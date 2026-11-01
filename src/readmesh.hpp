#include "classconfig.hpp"
#include "config.hpp"
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <string>

// 读取非结构网格
bool readmesh(const char* meshpath);
// 返回一个标识
void typetoken(char* a,char* b);

inline bool readmesh(const char* meshpath){
    printf("the path of mesh: %s\n", meshpath);
    FILE* meshfp; 
    int err = fopen_s(&meshfp, meshpath, "r");
    if (meshfp == NULL || err != 0 ){printf("not found mesh!\n"); return false;}
    if(fscanf_s(meshfp,"%d%d%d%d",&cc::node_num,&cc::face_num,&cc::cell_num,&cc::type_num)!=4)return false;
    printf("%d %d %d %d\n",cc::node_num,cc::face_num,cc::cell_num,cc::type_num);
    char buf[100];
    fgets(buf,sizeof(buf),meshfp);
    cc::NodeList.resize(cc::node_num);cc::FaceList.resize(cc::face_num);
    cc::CellList.resize(cc::cell_num);cc::bounds.reserve(cc::type_num-1);
    for(int i=1;i<=cc::type_num;i++){
        char name_before[30] = {}, name_after[30] = {}, dyh[2] = {};
        fscanf_s(meshfp, "%29s %1s %29s",
                 name_before, (unsigned)sizeof(name_before),
                 dyh, (unsigned)sizeof(dyh),
                 name_after, (unsigned)sizeof(name_after));
        printf("%s = %s\n",name_before,name_after);
        typetoken(name_before, name_after);
    }
    
    fgets(buf, sizeof(buf), meshfp);
    fgets(buf, sizeof(buf), meshfp);
    for(int i=0;i<cc::node_num;i++){
        int index;double x,y,z;
        fscanf_s(meshfp,"%d %lf %lf %lf",&index,&x,&y,&z);
        cc::NodeList[i] = cc::node_class(index,x,y,z);
    }
    fgets(buf, sizeof(buf), meshfp);
    fgets(buf, sizeof(buf), meshfp);
    printf("%s",buf);
    fclose(meshfp);
    return true;
}

inline void typetoken(char *a, char *b){
    if(strcmp(b,"WALL") == 0){cc::bounds.emplace_back(std::string(a),cc::WALL);}
    else if(strcmp(b,"FAR") == 0){cc::bounds.emplace_back(std::string(a),cc::FAR);}
    else if(strcmp(b,"SYM") == 0){cc::bounds.emplace_back(std::string(a),cc::SYM);}
}