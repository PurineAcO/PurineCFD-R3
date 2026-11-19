#pragma once

#include "data_structure/Vec5.h"

namespace JST {

struct dissipation {
    double Y = 0.0;
    cc::vec5 L;
    cc::vec5 Fd;
};

}
