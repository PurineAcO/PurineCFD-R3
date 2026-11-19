#pragma once
#include <cmath>

namespace cc {

struct vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
    vec3() = default;
    vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    vec3& operator+=(const vec3& o) {x += o.x;y += o.y;z += o.z;return *this;}
    vec3& operator-=(const vec3& o) {x -= o.x;y -= o.y;z -= o.z;return *this;}
    vec3& operator*=(double s) {x *= s;y *= s;z *= s;return *this;}
    void clear();
    double norm() const;
};

vec3 operator+(vec3 a, const vec3& b);
vec3 operator-(vec3 a, const vec3& b);
vec3 operator*(vec3 a, double s);
vec3 operator*(double s, vec3 a);
// 单位化
vec3 unit(const vec3& a);
// 点乘
double dot(const vec3& a, const vec3& b);
// 叉乘
vec3 cross(const vec3& a, const vec3& b);

inline void vec3::clear() {x = 0.0;y = 0.0; z = 0.0;}

inline double vec3::norm() const {return std::sqrt(x * x + y * y + z * z);}

}
