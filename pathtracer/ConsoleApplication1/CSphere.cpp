#include "stdafx.h"
#include "CSphere.h"
#include "SRay.h"
#include "SHitRecord.h"

//-------------------------------------------------------------------------------
// Constructors
//-------------------------------------------------------------------------------
CSphere::CSphere(Eigen::Vector3f center, float radius)
    : _c(center), _r(radius), _color(Eigen::Vector4f(1.0f, 0.0f, 0.0f, 0.0f))
{}

CSphere::CSphere(Eigen::Vector3f center, float radius, Eigen::Vector4f color)
    : _c(center), _r(radius), _color(color)
{}
//-------------------------------------------------------------------------------
// Intersection with a ray
//-------------------------------------------------------------------------------
void CSphere::intersect(const SRay& ray, SHitRecord& hit) const
{
    Eigen::Vector3f ec = ray.pos - _c;
    float A = ray.dir.dot(ray.dir);
    float B = (2 * ray.dir).dot(ec);
    float C = ec.dot(ec) - _r * _r;

    float discriminant = B * B - 4.0f * A * C;
    if (discriminant <= 0.001f)
        return;
    
    float sqrtDiscriminant = std::sqrtf(discriminant);
    float localT0 = (-B + sqrtDiscriminant) / (2 * A);
    float localT1 = (-B - sqrtDiscriminant) / (2 * A);
    float t = std::min(localT0, localT1);
    // fill struct
    hit.t1 = std::min<float>(hit.t1, t);
    hit.hittedObject = true;
}

//-------------------------------------------------------------------------------
// Normal at point P
//-------------------------------------------------------------------------------
Eigen::Vector3f CSphere::normalAt(const Eigen::Vector3f& p) const
{
    return (p - _c).normalized();
}

//-------------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------------
CSphere::~CSphere()
{}
