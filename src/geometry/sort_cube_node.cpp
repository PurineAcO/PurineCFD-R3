#include "geometry.h"


/*
为了方便计算多面体的体积,必须将各个点进行排序,使之满足某种关系
*/

std::vector<cc::node_class*> sort_cube_node(cc::cell_class &c){
    const cc::face_class* bottom = c.faces[0];
    const cc::face_class* top = c.faces[1];
    int shift = 0, flip = 0;
    double best = -1.0;
    for(int f = 0; f < 2; f++){
        for(int s = 0; s < 4; s++){
            double cost = 0.0;
            for(int i = 0; i < 4; i++){
                const int j = f ? (s - i + 8) % 4 : (s + i) % 4;
                const cc::vec3 d = top->node[j]->p - bottom->node[i]->p;
                cost += dot(d, d);
            }
            if(best < 0.0 || cost < best){best = cost; shift = s; flip = f;}
        }
    }
    std::vector<cc::node_class*> node_sorted(8, nullptr);
    for(int i = 0; i < 4; i++){
        node_sorted[i] = bottom->node[i];
        node_sorted[4 + i] = top->node[flip ? (shift - i + 8) % 4 : (shift + i) % 4];
    }
    return node_sorted;
}