#include "cell_class.h"
#include "face_class.h"

namespace cc {

cell_class::cell_class(int index_, int f1_, int f2_, int f3_, int f4_,
                       int f5_, int f6_, int fcnt_):
    index(index_), fcnt(fcnt_){
        faces[0] = &FaceList[f1_-1];faces[1] = &FaceList[f2_-1];
        faces[2] = &FaceList[f3_-1];faces[3] = &FaceList[f4_-1];
        faces[4] = &FaceList[f5_-1];faces[5] = &FaceList[f6_-1];
}

}
