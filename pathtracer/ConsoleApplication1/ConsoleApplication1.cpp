#include "stdafx.h"
#include <iostream>
#include <memory>
//#include <Eigen/Dense>
#include "CPNGImage.h"
#include "CSphere.h"
#include "CCamera.h"
#include "SRay.h"
#include "SHitRecord.h"


Eigen::Vector4f lerp(Eigen::Vector4f A, Eigen::Vector4f B, float t)
{
    return (1.0f - t) * A + t * B;
}


int main()
{
    Eigen::Vector4f red(1.0f, 0.0f, 0.0f, 1.0f);
    Eigen::Vector4f black(0.0f, 0.0f, 0.0f, 1.0f);
    Eigen::Vector4f blue(0.0f, 0.0f, 1.0f, 1.0f);

    CCamera camera(Eigen::Vector3f(0.0f, 0.0f, 5.0f), Eigen::Vector3f(0.0f, 0.0f, 0.0f));
	//std::unique_ptr<Surface> surface;
    CSphere sphere(Eigen::Vector3f(0.0f, 0.0f, -1.0f), 5.f, red);

    Eigen::Vector4f up(0.0f, 1.0f, 0.0f, 0.0f);
	Eigen::Vector4f right(1.0f, 0.0f, 0.0f, 0.0f);

	//generate some image
	const unsigned int width = 512, height = 512;
	CPNGImage img("Ball_2.png", width, height);
    int yx = width * height;
	for (unsigned y = 0; y < height; y++)
		for (unsigned x = 0; x < width; x++)
		{
            SRay primaryRay = camera.rayTrace(x, y, width, height);
            SHitRecord hit;
            sphere.intersect(primaryRay, hit);
            //std::cout << primaryRay.dir[0] << primaryRay.dir[1]<< primaryRay.dir[2] << std::endl;
            if (hit.hittedObject)
            {
                img.writeVectorToPixel(red,x, y);
            }
            else
            {
                Eigen::Vector4f background = lerp(black, blue, y/512.0f);
                img.writeVectorToPixel(background, x, y);
            }
		}
	img.encodeOneStep();

	std::cout << up << std::endl;
	std::cout << right << std::endl;
	std::cout << "dot:" << up.dot(right) << std::endl;
	std::cout << "back:" <<up.cross3(right) << std::endl;

	std::cout << "Epsilon:" << CSurface::Epsilon() << std::endl;
    return 0;
}

