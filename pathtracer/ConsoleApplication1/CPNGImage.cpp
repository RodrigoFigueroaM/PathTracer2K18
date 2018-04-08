#include "stdafx.h"
#include "CPNGImage.h"

//-------------------------------------------------------------------------------
// Constructors
//-------------------------------------------------------------------------------
CPNGImage::CPNGImage()
	: CPNGImage("default.png", 100, 100) {}

CPNGImage::CPNGImage(const std::string name, const unsigned int & width, const unsigned int & height)
	: _name(name), _width(width), _height(height) 
{
	const unsigned int bytesPerPixel = 4;
	_image.resize(width * height * bytesPerPixel);
}
//-------------------------------------------------------------------------------
// Write Eigen vector to pixel with x,y coordinates 
//-------------------------------------------------------------------------------
void  CPNGImage::writeVectorToPixel(const Eigen::Vector4f & V, const unsigned int x, const unsigned y)
{
	unsigned int pixel = 4 * _width * y + 4 * x;
	_image[pixel + 0] = V.x() * 255;
	_image[pixel + 1] = V.y() * 255;
	_image[pixel + 2] = V.z() * 255;
	_image[pixel + 3] = 255;
}

//-------------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------------
CPNGImage::~CPNGImage() {}

//-------------------------------------------------------------------------------
// Private members
//-------------------------------------------------------------------------------
void CPNGImage::encodeOneStep()
{
	//Encode the image
	unsigned error = lodepng::encode(_name, _image, _width, _height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}