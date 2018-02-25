// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>

int main()
{
	Eigen::Vector3f up(0.0, 1.0, 0.0);
	Eigen::Vector3f right(1.0, 0.0, 0.0);

	std::cout << up << std::endl;
	std::cout << right << std::endl;
	std::cout << "dot:" << up.dot(right) << std::endl;
	std::cout << "back:" << up.cross(right) << std::endl;
    return 0;
}

