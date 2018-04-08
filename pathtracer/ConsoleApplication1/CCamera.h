#pragma once
#include <Eigen/Geometry>

struct SRay;

class CCamera {
public:
    CCamera(Eigen::Vector3f position, Eigen::Vector3f lookAt, Eigen::Vector3f up = Eigen::Vector3f(0.0f, 1.0f, 0.0f));
    Eigen::Vector3f position() const;
    Eigen::Vector3f direction() const;
    Eigen::Vector3f up() const;
    SRay rayTrace(const int& x, const int& y, const int& width, const int& height, const float s = 1.0) const;
    ~CCamera();
private:
    Eigen::Vector3f _pos;
    Eigen::Vector3f _dir;
    Eigen::Vector3f _up;
    Eigen::Vector3f _right;
    Eigen::Vector3f _back;
};

inline Eigen::Vector3f CCamera::position() const {
    return _pos;
}


inline Eigen::Vector3f CCamera::direction() const {
    return _dir;
}


inline Eigen::Vector3f CCamera::up() const {
    return _up;
}