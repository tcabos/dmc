#ifndef WALKER
#define WALKER

#include "vec3.hpp"

struct alignas(64) Walker
{
    Vec3 e1;
    Vec3 e2;
    Walker(const Vec3 &e1_ = Vec3(), const Vec3 &e2_ = Vec3()) : e1(e1_), e2(e2_) {}
};

#endif 