
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SHAPES_CURVE_H
#define PBRT_SHAPES_CURVE_H

// shapes/curve.h*
#include "shape.h"

namespace pbrt {
struct CurveCommon;

// CurveType Declarations
enum class CurveType { Flat, Cylinder, Ribbon };

// CurveCommon Declarations
struct CurveCommon {
    CurveCommon(const Point3f c[4], Float w0, Float w1, CurveType type,
                const Normal3f* norm);
    const CurveType type;
    Point3f         cpObj[4];
    Float           width[2];
    Normal3f        n[2];
    Float           normalAngle, invSinNormalAngle;
};

// Curve Declarations
class Curve : public Shape {
public:
    // Curve Public Methods
    Curve(const Transform* ObjectToWorld, const Transform* WorldToObject,
          bool reverseOrientation, const std::shared_ptr<CurveCommon>& common,
          Float uMin, Float uMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
        common(common),
        uMin(uMin),
        uMax(uMax) {}
    Bounds3f ObjectBound() const;
    bool     Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
                       bool testAlphaTexture) const;

    bool IntersectAABB(const Ray& ray, Float* tHit, SurfaceInteraction* isect) const;
    bool IntersectPhantom(const Ray& ray, Float* tHit, SurfaceInteraction* isect) const;

    Float       Area() const;
    Interaction Sample(const Point2f& u, Float* pdf) const;

private:
    // Curve Private Methods
    bool recursiveIntersect(const Ray& r, Float* tHit,
                            SurfaceInteraction* isect, const Point3f cp[4],
                            const Transform& rayToObject, Float u0, Float u1,
                            int depth) const;

    // Curve Private Data
    const std::shared_ptr<CurveCommon> common;
    const Float uMin, uMax;
};

std::vector<std::shared_ptr<Shape>> CreateCurveShape(const Transform* o2w,
                                                     const Transform* w2o,
                                                     bool             reverseOrientation,
                                                     const ParamSet&  params);

/****************************************************************************************************/
// CurveCommon Declarations
struct QuadraticCurveCommon {
    QuadraticCurveCommon(const Point3f c[3], Float w0, Float w1, CurveType type,
                         const Normal3f* norm);
    const CurveType type;
    Point3f         cpObj[3];
    Float           width[2];
    Normal3f        n[2];
    Float           normalAngle, invSinNormalAngle;
};

// Curve Declarations
class QuadraticCurve : public Shape {
public:
    // Curve Public Methods
    QuadraticCurve(const Transform* ObjectToWorld, const Transform* WorldToObject,
                   bool reverseOrientation, const std::shared_ptr<QuadraticCurveCommon>& common,
                   Float uMin, Float uMax)
        : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
        common(common),
        uMin(uMin),
        uMax(uMax) {}
    Bounds3f ObjectBound() const;
    bool     Intersect(const Ray& ray, Float* tHit, SurfaceInteraction* isect,
                       bool testAlphaTexture) const;
    Float       Area() const;
    Interaction Sample(const Point2f& u, Float* pdf) const;

private:
    // Curve Private Methods
    bool recursiveIntersect(const Ray& r, Float* tHit,
                            SurfaceInteraction* isect, const Point3f cp[3],
                            const Transform& rayToObject, Float u0, Float u1,
                            int depth) const;

    // Curve Private Data
    const std::shared_ptr<QuadraticCurveCommon> common;
    const Float uMin, uMax;
};

/****************************************************************************************************/

struct RayConeIntersection {
    bool intersect(float r, float dr) {
        float r2  = r * r;
        float drr = r * dr;

        float ddd = cd.x * cd.x + cd.y * cd.y;
        if(ddd == 0) { ddd = 1; }

        dp = c0.x * c0.x + c0.y * c0.y;
        float cdd = c0.x * cd.x + c0.y * cd.y;
        float cxd = c0.x * cd.y - c0.y * cd.x;

        float c    = ddd;
        float b    = cd.z * (drr - cdd);
        float cdz2 = cd.z * cd.z;
        ddd += cdz2;
        float a = 2 * drr * cdd + cxd * cxd - ddd * r2 + dp * cdz2;

 #ifdef KEEP_DR2
        float qs = (dr * dr) / ddd;
        a -= qs * cdd * cdd;
        b -= qs * cd.z * cdd;
        c -= qs * cdz2;
 #endif

        // We will add c0 . z to s and splatter if needed
        float det = b * b - a * c;
        s   = (b - (det > 0 ? sqrt(det) : 0)) / c;
        dt  = (s * cd.z - cdd) / ddd;
        dc  = s * s + dp;
        sp  = cdd / cd.z;
        dp += sp * sp;

        return det > 0;
    }

    Vector3f c0;
    Vector3f cd;
    float    s;
    float    dt;
    float    dp;
    float    dc;
    float    sp;
};
} // namespace pbrt

#endif // PBRT_SHAPES_CURVE_H
