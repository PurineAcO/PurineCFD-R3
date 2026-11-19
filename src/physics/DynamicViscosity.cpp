#include "physics/Physics.h"

double ph::dynamic_viscosity(double T) {
    return mu0 * (T / T0) * std::sqrt(T / T0) * (T0 + Ts) / (T + Ts);
}
