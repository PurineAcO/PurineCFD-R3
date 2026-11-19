#include "physics/Physics.h"

double ph::_MUT(double mu, double T) {
    return mu * (T + 3 * Ts) / (2 * T * (T + Ts));
}
