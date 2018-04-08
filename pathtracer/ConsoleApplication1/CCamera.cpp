#include "stdafx.h"
#include "CCamera.h"
#include "SRay.h"

CCamera::CCamera(Eigen::Vector3f position,
        Eigen::Vector3f lookAt, 
        Eigen::Vector3f up) 
    : _pos(position)
{
    // uvw camera local space
    //Eigen::Vector3f tempPos([0], position[1], position[2]);
    _dir = (lookAt - _pos).normalized();
    _back = (_pos - lookAt).normalized();
    _right = up.cross(_dir);
    _up = _back.cross(_right);
    //_center = _pos + _dir;
}

SRay CCamera::rayTrace(const int& x, const int& y, const int& width, const int& height, const float s) const
{
    //Eigen::Vector3f tempPos(_pos[0], _pos[1], _pos[2]);
    int imageY = height - y;
    float sampleX = s * (x - 0.5f * (width + 1.0f));
    float sampleY = s * (imageY - 0.5f * (height + 1.0f));
    float camDistanceToViewingPlane = 1.0f;
    Eigen::Vector3f sw = sampleX * _right +  sampleY * _up - _back * camDistanceToViewingPlane;
    Eigen::Vector3f swnorm = (sw - _pos).normalized();
    SRay camray;
    camray.pos = _pos;
    camray.dir = swnorm;
    //Eigen::Vector4f(swnorm[0], swnorm[1], swnorm[2], 0.0f);
    return camray;
}

CCamera::~CCamera()
{
}
