// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>
#include "lodepng.h"

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

void writeVectorToPixel(const Eigen::Vector3f & V, std::vector<unsigned char> & image, const unsigned int pixel)
{
	image[pixel + 0] = V.x() * 255;
	image[pixel + 1] = V.y() * 255;
	image[pixel + 2] = V.z() * 255;
	image[pixel + 3] = 255;
}

int main()
{


	Eigen::Vector3f up(0.0, 1.0, 0.0);
	Eigen::Vector3f right(1.0, 0.0, 0.0);

	const char* filename = "test.png";

	//generate some image
	const unsigned int width = 512, height = 512;
	std::vector<unsigned char> image;
	const unsigned int bytesPerPixel = 4;
	image.resize(width * height * bytesPerPixel );
	for (unsigned y = 0; y < height; y++)
		for (unsigned x = 0; x < width; x++)
		{
			unsigned int pixel = bytesPerPixel * width * y + bytesPerPixel * x;
			writeVectorToPixel(right, image, pixel);
		}
	encodeOneStep(filename, image, width, height);

	std::cout << up << std::endl;
	std::cout << right << std::endl;
	std::cout << "dot:" << up.dot(right) << std::endl;
	std::cout << "back:" << up.cross(right) << std::endl;

    return 0;
}

