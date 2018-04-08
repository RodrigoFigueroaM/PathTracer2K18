#pragma once

#include <Eigen/Dense>

struct SRay;
struct SHitRecord;

class CSurface
{
public:
    static const float& Epsilon();
    virtual void intersect(const SRay& ray, SHitRecord& hit) const = 0;
    virtual Eigen::Vector3f normalAt(const Eigen::Vector3f& V) const = 0;
    virtual Eigen::Vector4f color() const;
};

//-------------------------------------------------------------------------------
// static members
//-------------------------------------------------------------------------------
inline const float& CSurface::Epsilon()
{
	return 0.001f;
}


inline Eigen::Vector4f CSurface::color() const
{
    return Eigen::Vector4f (1.0f, 1.0f, 1.0f, 1.0f);
}
