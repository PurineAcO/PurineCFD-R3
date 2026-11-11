#include "face_class.h"
#include "cell_class.h"
#include "node_class.h"

namespace cc {

face_class::face_class(int index_, int p1_, int p2_, int p3_, int p4_,
                       int c1_, int c2_, short type_, short pcnt_):
    index(index_), type(type_), pcnt(pcnt_){
        node[0] = &NodeList[p1_-1];node[1] = &NodeList[p2_-1];
        node[2] = &NodeList[p3_-1];node[3] = &NodeList[p4_-1];
        nei[0] = c1_ > 0 ? &CellList[c1_-1] : nullptr;
        nei[1] = c2_ > 0 ? &CellList[c2_-1] : nullptr;
}

}
