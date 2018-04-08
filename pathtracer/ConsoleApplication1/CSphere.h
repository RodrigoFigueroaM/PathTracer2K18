#pragma once
#include "CSurface.h"

class CSphere : public CSurface
{
public:
    CSphere(Eigen::Vector3f center, float radius);
    CSphere(Eigen::Vector3f center, float radius, Eigen::Vector4f color);
	Eigen::Vector3f center() const;
	float radius() const;
    void intersect(const SRay& ray, SHitRecord& hit) const override;
    Eigen::Vector3f normalAt(const Eigen::Vector3f& P) const override;
    Eigen::Vector4f color() const override;
	~CSphere();

private:
    Eigen::Vector4f _color;
	Eigen::Vector3f _c;
	float _r;

};

//-------------------------------------------------------------------------------
// accessors
//-------------------------------------------------------------------------------
inline Eigen::Vector3f CSphere::center() const 
{
    return _c;
}

inline float CSphere::radius() const
{
    return _r;
}

inline Eigen::Vector4f CSphere::color() const
{
    return _color;
}
