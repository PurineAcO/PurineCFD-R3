#pragma once

namespace cc {

struct vec5 {
    double c = 0.0, x = 0.0, y = 0.0, z = 0.0, e = 0.0;
    vec5() = default;
    vec5(double c_, double x_, double y_, double z_, double e_): c(c_), x(x_), y(y_), z(z_), e(e_) {}
    vec5& operator+=(const vec5& o) {c += o.c;x += o.x;y += o.y;z += o.z;e += o.e;return *this;}
    vec5& operator-=(const vec5& o) {c -= o.c;x -= o.x;y -= o.y;z -= o.z;e -= o.e;return *this;}
    vec5& operator*=(double k) {c *= k;x *= k; y *= k;z *= k;e *= k;return *this;}
    double& operator[](int i);
    double operator[](int i) const;
    void clear();
};

vec5 operator+(vec5 a, const vec5& b);
vec5 operator-(vec5 a, const vec5& b);
vec5 operator*(vec5 a, double s);
vec5 operator*(double s, vec5 a);

inline void vec5::clear() {c = 0.0;x = 0.0;y = 0.0;z = 0.0;e = 0.0;}

}
