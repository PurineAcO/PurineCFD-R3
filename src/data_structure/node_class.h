#pragma once
#include "classconfig.h"

namespace cc {
// 节点
struct node_class {
    int index = 0;    // 节点编号
    vec3 p;           // 节点坐标             
    node_class() = default;
    node_class(int number_, double x_, double y_, double z_);
};
inline node_class::node_class(int number_,double x_,double y_,double z_):index(number_),p(vec3(x_,y_,z_)){}
} 