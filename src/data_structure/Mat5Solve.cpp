#include <cmath>

#include "data_structure/Mat5.h"

namespace cc {

bool mat5::solve(vec5& b) const {
    double m[5][6];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) m[i][j] = a[i][j];
        m[i][5] = (&b.c)[i];
    }
    for (int k = 0; k < 5; k++) {
        int p = k;
        for (int i = k + 1; i < 5; i++) {
            if (std::abs(m[i][k]) > std::abs(m[p][k])) p = i;
        }
        if (p != k) {
            for (int j = 0; j < 6; j++) {
                const double t = m[k][j];
                m[k][j] = m[p][j];
                m[p][j] = t;
            }
        }
        if (std::abs(m[k][k]) < 1e-14) m[k][k] = m[k][k] >= 0.0 ? 1e-14 : -1e-14;
        for (int i = k + 1; i < 5; i++) {
            const double f = m[i][k] / m[k][k];
            for (int j = k; j < 6; j++) m[i][j] -= f * m[k][j];
        }
    }
    double x[5];
    for (int i = 4; i >= 0; i--) {
        double s = m[i][5];
        for (int j = i + 1; j < 5; j++) s -= m[i][j] * x[j];
        x[i] = s / m[i][i];
    }
    b = vec5(x[0], x[1], x[2], x[3], x[4]);
    return true;
}

}
