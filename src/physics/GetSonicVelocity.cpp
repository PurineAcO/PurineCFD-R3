#include "physics/Physics.h"

double ph::get_sonic_velocity(double T) {
    return std::sqrt(ph::gamma * ph::R * T);
}
