#pragma once
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include "lodepng.h"


//TODO: forward declare Eigen

class CPNGImage
{
public:
	CPNGImage();
	CPNGImage(const std::string name, const unsigned int & width, const unsigned int & height);
	const std::string name();
	const unsigned int width();
	const unsigned int height();
	void writeVectorToPixel(const Eigen::Vector4f & V, const unsigned int x, const unsigned y);
	void encodeOneStep();
	~CPNGImage();
private:
	
	std::vector<unsigned char> _image;
	const std::string _name;
	const unsigned int _width;
	const unsigned int _height;
};

//-------------------------------------------------------------------------------
// accessors
//-------------------------------------------------------------------------------
inline const std::string CPNGImage::name()
{
	return _name;
}

inline const unsigned int CPNGImage::width()
{
	return _width;
}

inline const unsigned int CPNGImage::height()
{
	return _height;
}
