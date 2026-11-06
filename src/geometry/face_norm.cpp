#include "geometry.h"

/*
return normal vector of face 
Warning:
1. 在执行这个函数之前,请务必保证四边形(如果提供的是)四个点是共面的,如果发生了微微的翘曲,请执行修复函数
2. 返回的面法向是带面积的单位向量
*/

cc::vec3 face_norm(cc::face_class &f){
    if(f.pcnt == 3) return 0.5*cc::cross(f.node[1]->p-f.node[0]->p, f.node[2]->p-f.node[0]->p);
    return 0.5*cc::cross(f.node[2]->p-f.node[0]->p, f.node[3]->p-f.node[1]->p);
}