#include "cell_class.h"
#include "face_class.h"

namespace cc {

cell_class::cell_class(int index_, int f1_, int f2_, int f3_, int f4_,
                       int f5_, int f6_, int fcnt_):
    index(index_), fcnt(fcnt_){
        const int ids[6] = {f1_, f2_, f3_, f4_, f5_, f6_};
        for(int j = 0; j < fcnt_ && j < 6; j++){
            faces[j] = ids[j] > 0 ? &FaceList[ids[j]-1] : nullptr;
        }
}

}
