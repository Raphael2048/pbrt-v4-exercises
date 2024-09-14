// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <pbrt/util/math.h>
#include <pbrt/util/print.h>
#include <pbrt/util/stats.h>
#include <pbrt/util/transform.h>
#include <pbrt/util/vecmath.h>
#include <pbrt/util/check.h>

#include <algorithm>
#include <cmath>
#include <type_traits>

namespace pbrt {

template <>
std::string internal::ToString2<Interval>(Interval x, Interval y) {
    return StringPrintf("[ %s %s ]", x, y);
}

template <>
std::string internal::ToString3<Interval>(Interval x, Interval y, Interval z) {
    return StringPrintf("[ %s %s %s ]", x, y, z);
}

template <typename T>
std::string internal::ToString2(T x, T y) {
    if (std::is_floating_point_v<T>)
        return StringPrintf("[ %f, %f ]", x, y);
    else
        return StringPrintf("[ %d, %d ]", x, y);
}

template <typename T>
std::string internal::ToString3(T x, T y, T z) {
    if (std::is_floating_point_v<T>)
        return StringPrintf("[ %f, %f, %f ]", x, y, z);
    else
        return StringPrintf("[ %d, %d, %d ]", x, y, z);
}

template std::string internal::ToString2(float, float);
template std::string internal::ToString2(double, double);
template std::string internal::ToString2(int, int);
template std::string internal::ToString3(float, float, float);
template std::string internal::ToString3(double, double, double);
template std::string internal::ToString3(int, int, int);

// Quaternion Method Definitions
std::string Quaternion::ToString() const {
    return StringPrintf("[ %f, %f, %f, %f ]", v.x, v.y, v.z, w);
}

constexpr bool USE_PRECISE_CHECK = true;

// Points from a cube vertices can't be same line, and must can be solved
Point3f SolveIntersectPoint(Float* A, Float* B, Float* C) {
    if (A[0] == 0)
    {
        if (B[0] != 0.0) std::swap(A, B);
        if (C[0] != 0.0) std::swap(A, C);
    }

    Point3f Result;
    Float k1 = B[0] / A[0];
    Float k2 = C[0] / A[0];
    for (int i = 0; i < 4; ++i) {
        B[i] -= k1 * A[i];
        C[i] -= k2 * A[i];
    }

    if (B[1] == 0) {
        std::swap(B, C);
    }
    Float k = C[1] / B[1];
    for (int i = 1; i < 4; ++i) {
        C[i] -= B[i] * k;
    }

    C[3] /= C[2];
    C[2] = 1.0;
    Result.z = C[3];

    B[3] -= Result.z * B[2];
    B[2] = 0;
    B[3] /= B[1];
    B[1] = 1.0f;
    Result.y = B[3];

    A[3] -= Result.z * A[2];
    A[2] = 0;
    A[3] -= Result.y * A[1];
    A[1] = 0;
    A[3] /= A[0];
    A[0] = 1.0f;
    Result.x = A[3];
    return Result;
}

// 给出三个点，得到外接圆
Point3f CircumCircle(const Point3f& a, const Point3f& b, const Point3f& c)
{
    Vector3f t1 = b - a;
    Float HalfPlaneA[4];
    HalfPlaneA[0] = t1.x;
    HalfPlaneA[1] = t1.y;
    HalfPlaneA[2] = t1.z;
    HalfPlaneA[3] = Dot(Vector3f(a + b) * 0.5f, t1);

    Vector3f t2 = c - a;

    Float HalfPlaneB[4];
    HalfPlaneB[0] = t2.x;
    HalfPlaneB[1] = t2.y;
    HalfPlaneB[2] = t2.z;
    HalfPlaneB[3] = Dot(Vector3f(a + c) * 0.5f, t2);

    Vector3f n = Cross(t1, t2);

    Float TriPlane[4];
    TriPlane[0] = n.x;
    TriPlane[1] = n.y;
    TriPlane[2] = n.z;
    TriPlane[3] = Dot(Vector3f(a), n);

    Point3 IntersectP = SolveIntersectPoint(HalfPlaneA, HalfPlaneB, TriPlane);

    // float d1 = DistanceSquared(a, IntersectP);
    // float d2 = DistanceSquared(b, IntersectP);
    // float d3 = DistanceSquared(c, IntersectP);
    // DCHECK(d1 - d2 < 1e-3);
    // DCHECK(d2 - d3 < 1e-3);

    return IntersectP;
}

DirectionCone BoundSubtendedDirections(const Bounds3f& b, Point3f p) {
    if constexpr (USE_PRECISE_CHECK) {
        if (Inside(p, b)) {
            return DirectionCone::EntireSphere();
        }
        else {
            Vector3f Directions[8];
            for (int i = 0; i < 8; ++i)
            {
                Directions[i] = Normalize(b.Corner(i) - p);
            }

            Vector3f w(0, 0, 1);
            Float cosTheta = 1;
            w = Directions[0];
            
            for (int i = 1; i < 8; ++i) {
                if (Dot(Directions[i], w) < cosTheta) {
                    w = Directions[i];
                    cosTheta = 1;
                    for (int j = 0; j < i; ++j) {
                        if (Dot(Directions[j], w) < cosTheta) {
                            w = Normalize(Directions[i] + Directions[j]);
                            cosTheta = Dot(w, Directions[i]);
                            for (int k = 0; k < j; ++k) {
                                if (Dot(Directions[k], w) < cosTheta) {
                                    Point3f center =
                                        CircumCircle(Point3f(Directions[i]),
                                                    Point3f(Directions[j]),
                                                      Point3f(Directions[k]));
                                    w = Normalize(Vector3f(center));
                                    cosTheta = Dot(w, Directions[i]);
                                }
                            }
                        }
                    }
                }
            }
            return DirectionCone(w, cosTheta);
        }
    } else {
        // Compute bounding sphere for _b_ and check if _p_ is inside
        Float radius;
        Point3f pCenter;
        b.BoundingSphere(&pCenter, &radius);
        if (DistanceSquared(p, pCenter) < Sqr(radius))
            return DirectionCone::EntireSphere();

        // Compute and return _DirectionCone_ for bounding sphere
        Vector3f w = Normalize(pCenter - p);
        Float sin2ThetaMax = Sqr(radius) / DistanceSquared(pCenter, p);
        Float cosThetaMax = SafeSqrt(1 - sin2ThetaMax);
        return DirectionCone(w, cosThetaMax);
    }
}


// DirectionCone Function Definitions
DirectionCone Union(const DirectionCone &a, const DirectionCone &b) {
    // Handle the cases where one or both cones are empty
    if (a.IsEmpty())
        return b;
    if (b.IsEmpty())
        return a;

    // Handle the cases where one cone is inside the other
    Float theta_a = SafeACos(a.cosTheta), theta_b = SafeACos(b.cosTheta);
    Float theta_d = AngleBetween(a.w, b.w);
    if (std::min(theta_d + theta_b, Pi) <= theta_a)
        return a;
    if (std::min(theta_d + theta_a, Pi) <= theta_b)
        return b;

    // Compute the spread angle of the merged cone, $\theta_o$
    Float theta_o = (theta_a + theta_d + theta_b) / 2;
    if (theta_o >= Pi)
        return DirectionCone::EntireSphere();

    // Find the merged cone's axis and return cone union
    Float theta_r = theta_o - theta_a;
    Vector3f wr = Cross(a.w, b.w);
    if (LengthSquared(wr) == 0)
        return DirectionCone::EntireSphere();
    Vector3f w = Rotate(Degrees(theta_r), wr)(a.w);
    return DirectionCone(w, std::cos(theta_o));
}

std::string DirectionCone::ToString() const {
    return StringPrintf("[ DirectionCone w: %s cosTheta: %f ]", w, cosTheta);
}

}  // namespace pbrt
