#pragma once
#include <memory>

class CSurface;
struct SHitRecord
{
    bool hittedObject = false;
    float t0;
    float t1;
    std::unique_ptr<CSurface> object;
};