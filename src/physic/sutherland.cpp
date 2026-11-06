#include "physic.h"

double ph::dynamic_viscosity(double T) {
    return mu0 * (T / T0) * std::sqrt(T / T0) * (T0 + Ts) / (T + Ts);
}

double ph::_MUT(double mu, double T) {
    return mu * (T + 3 * Ts) / (2 * T * (T + Ts));
}