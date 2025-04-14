#pragma once

#include <cmath>

class alignas(32) Vec3
{
public:
    double x, y, z;
    Vec3(double x_ = 0, double y_ = 0, double z_ = 0) : x(x_), y(y_), z(z_) {}
    inline double norm() const { return std::sqrt(x * x + y * y + z * z); }
    inline double distance(const Vec3 &other) const
    {
        double dx = x - other.x, dy = y - other.y, dz = z - other.z;
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
};
