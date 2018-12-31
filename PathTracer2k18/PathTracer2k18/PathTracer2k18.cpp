// PathTracer2k18.cpp 

#include "pch.h"
#include <stdio.h>
#include <ctime>
#include <limits>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

struct Vector3f {
	float x, y, z;
	Vector3f(float x0 = 0.f, float y0 = 0.f, float z0 = 0.f) { x = x0; y = y0; z = z0; }
	Vector3f operator+(const Vector3f &b) const { return Vector3f(x + b.x, y + b.y, z + b.z); }
	Vector3f operator-(const Vector3f &b) const { return Vector3f(x - b.x, y - b.y, z - b.z); }
	Vector3f operator*(float b) const { return Vector3f(x*b, y*b, z*b); }
	Vector3f operator/(float b) const { return Vector3f(x / b, y / b, z / b); }
	Vector3f Mult(const Vector3f &b) const { return Vector3f(x*b.x, y*b.y, z*b.z); }
	Vector3f& Norm() { return *this = *this * (1 / sqrt(x*x + y * y + z * z)); }
	float Length() { return sqrt(x*x + y * y + z * z); }
	float Dot(const Vector3f &b) const { return x * b.x + y * b.y + z * b.z; }
	Vector3f operator%(const Vector3f &b) const { return Vector3f(y*b.z - z * b.y, z*b.x - x * b.z, x*b.y - y * b.x); }
	const float& operator[](size_t i) const { return i == 0 ? x : (i == 1 ? y : z); }
};

std::mt19937 mersenneTwister;
std::uniform_real_distribution<float> uniform;

#define RND (2.0f*uniform(mersenneTwister)-1.0f)
#define RND2 (uniform(mersenneTwister))

namespace Globals {
	const float epsilon = 1e-6f;
	const int width = 512;
	const int height = 512;
	const float pi = 3.1415926535f;
	const int maxdepth = 5;
	const char* filename = "PathTracer2k18.ppm";
	const float spp = 750.f;
}

struct Ray {
	Vector3f origin;
	Vector3f dir;
};

class Hitable {
public:
	enum class Type : int {
		Diffuse,
		Specular, // reflective
		Refraction,
	};

	// -- Data
	Vector3f color;
	float emission;
	Type type;

	// -- Interface
	virtual float Intersect(const Ray& ray) const = 0;
	virtual Vector3f Normal(const Vector3f& point) const = 0;

	// -- Constructors
	Hitable() : color(Vector3f(255.f, 0.f, 255.f)), emission(0.0f), type(Type::Diffuse) {
	}

	Hitable(Vector3f color, float emission, Type type) :
		color(color), emission(emission), type(type) {
	}

	~Hitable() = default;
};

class Sphere : public Hitable {
public:
	Vector3f center;
	float radius;

	Sphere(float radius, Vector3f center) :
		radius(radius), center(center) {
	}
	Sphere(float radius, Vector3f center, Vector3f color, float emission, Hitable::Type type) :
		radius(radius), center(center), Hitable(color, emission, type) {
	}
	~Sphere() = default;

	float Intersect(const Ray& ray) const {
		float b = ((ray.origin - center) * 2.0f).Dot(ray.dir);
		float c = (ray.origin - center).Dot((ray.origin - center)) - (radius * radius);
		float discriminant = b * b - (4.0f*c);
		if (discriminant < 0.0f)
			return 0.0f;
		else
			discriminant = sqrt(discriminant);

		float x1 = -b + discriminant;
		float x2 = -b - discriminant;

		return (x2 > Globals::epsilon) ? x2 / 2.0f : (x1 > Globals::epsilon) ? x1 / 2.0f : 0.0f;
	};

	Vector3f Normal(const Vector3f& point) const {
		return (point - center).Norm();
	}
};

class Plane : public Hitable {
public:
	Vector3f normal;
	float distance;
	Plane(float distance, Vector3f normal) :
		distance(distance), normal(normal) {
	}
	Plane(float distance, Vector3f normal, Vector3f color, float emission, Hitable::Type type) :
		distance(distance), normal(normal), Hitable(color, emission, type) {
	}

	float Intersect(const Ray& ray) const {
		float d = normal.Dot(ray.dir);
		if (d != 0.0f) {
			float t = -1.0f * ((normal.Dot(ray.origin) + distance) / d);
			return (t > Globals::epsilon) ? t : 0.0f;
		}
		return 0.0f;
	}
	Vector3f Normal(const Vector3f& point) const {
		return normal;
	}
};

struct Intersection {
	float t = FLT_MAX;
	Hitable*  object = nullptr;
};

class Scene {
private:
	std::vector<Hitable*> objects;

public:
	void Add(Hitable* object) {
		objects.emplace_back(object);
	}

	Intersection Intersect(const Ray& ray) const {
		Intersection closestintersection;

		for (auto object : objects) {
			float t = object->Intersect(ray);
			if (t > Globals::epsilon && t < closestintersection.t) {
				closestintersection.t = t;
				closestintersection.object = object;
			}
		}

		return closestintersection;
	}
};

Vector3f CamCords(const float x, const float y) {
	float w = Globals::width;
	float h = Globals::height;
	float fovx = Globals::pi / 4;
	float fovy = (h / w) * fovx;
	return Vector3f(((2.f*x - w) / w) * tan(fovx),
		-((2.f*y - h) / h) * tan(fovy),
		-1.0);
}

// create orthonormal system
void OrthonormalSystem(Vector3f& a, Vector3f& b, Vector3f& c) {
	if (std::abs(a.x) > std::abs(a.y)) {
		// project to the y = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.0f / sqrt(a.x* a.x + a.z * a.z);
		b = Vector3f(-a.z * invLen, 0.0f, a.x * invLen);
	}
	else {
		// project to the x = 0 plane and construct a normalized orthogonal vector in this plane
		float invLen = 1.0f / sqrt(a.y * a.y + a.z * a.z);
		b = Vector3f(0.0f, a.z * invLen, -a.y * invLen);
	}
	c = Vector3f(a.y * b.z - a.z * b.y,
		a.z* b.x - a.x * b.z,
		a.x * b.y - a.y * b.x); // cross
}

// Class for generating the Halton low-discrepancy series for Quasi
// Monte Carlo integration.
class Halton {
	float value, inv_base;
public:
	void Number(int i, int base) {
		float f = inv_base = 1.0f / base;
		value = 0.0f;
		while (i > 0) {
			value += f * (float)(i%base);
			i /= base;
			f *= inv_base;
		}
	}
	void Next() {
		float r = 1.0f - value - 0.0000001f;
		if (inv_base < r) value += inv_base;
		else {
			float h = inv_base, hh;
			do { hh = h; h *= inv_base; } while (h >= r);
			value += hh + h - 1.0f;
		}
	}
	float Get() { return value; }
};

// Uniform sampling on a hemisphere from
// http://www.rorydriscoll.com/2009/01/07/better-sampling/
Vector3f Hemisphere(const float u, const float v) {
	const float r = sqrtf(1.0f - u * u);
	const float phi = 2.0f * Globals::pi * v;

	return Vector3f(cos(phi) * r, sin(phi) * r, u);
}


Vector3f Reflect(const Vector3f& incident, const Vector3f& normal) {
	return (incident - normal * (incident.Dot(normal) * 2.0f)).Norm();
}

void Trace(Ray& ray, const Scene& scene, int depth, Vector3f& color, Halton& hal, Halton& hal2) {
	if (depth >= Globals::maxdepth)
		return;

	Intersection intersection = scene.Intersect(ray);
	if (!intersection.object)
		return;

	// get hit point and normal at that hit point
	Vector3f hitpoint = ray.origin + ray.dir * intersection.t;
	Vector3f hitpointnrml = intersection.object->Normal(hitpoint);

	// new ray but we can reuse like its bouncing ;)
	ray.origin = hitpoint;

	// compute emission L_e(x,w) part of the rendering equation 
	const float emission = intersection.object->emission;

	color = color + intersection.object->color * emission;

	// Diffuse BRDF
	auto& object = intersection.object;
	if (intersection.object->type == Hitable::Type::Diffuse) {
		Vector3f rotx, roty;
		OrthonormalSystem(hitpointnrml, rotx, roty);
		Vector3f sampleddir = Hemisphere(RND2, RND2);

		//Transform the Random Samples to the Shaded Point Local Coordinate System
		Vector3f rotateddir;
		rotateddir.x = Vector3f(rotx.x, roty.x, hitpointnrml.x).Dot(sampleddir);
		rotateddir.y = Vector3f(rotx.y, roty.y, hitpointnrml.y).Dot(sampleddir);
		rotateddir.z = Vector3f(rotx.z, roty.z, hitpointnrml.z).Dot(sampleddir);
		ray.dir = rotateddir;

		float cost = ray.dir.Dot(hitpointnrml);

		Vector3f tmpcolor;
		Trace(ray, scene, depth + 1, tmpcolor, hal, hal2);
		color = color + (tmpcolor.Mult(intersection.object->color)) * cost * 0.1f;
	}

	// Specular BRDF  the perfect reflection direction.
	if (intersection.object->type == Hitable::Type::Specular) {
		ray.dir = Reflect(ray.dir, hitpointnrml);
		Vector3f tmpcolor;
		Trace(ray, scene, depth + 1, tmpcolor, hal, hal2);
		color = color + tmpcolor;
	}

	// Glass/refractive BTDF  - we use the vector version of Snell's law and Fresnel's law
	// to compute the outgoing reflection and refraction directions and probability weights.
	if (intersection.object->type == Hitable::Type::Refraction) {
		float n = 1.5f;
		float r0 = (1.0f - n) / (n + 1.0f);
		r0 = r0 * r0;

		// ray coming from inside the medium
		if (hitpointnrml.Dot(ray.dir) > 0) {
			hitpointnrml = hitpointnrml * -1.0f;
			n = 1.0f / n;
		}

		n = 1.0f / n;
		float cost1 = hitpointnrml.Dot(ray.dir) * -1.0f;
		float cost2 = 1.0f - n * n *(1.0f - cost1 * cost1);
		auto SchlicksApproximation = [](const float& r0, const float& cost)->float {
			return r0 + (1.0f - r0) * ((1 - cost) * (1 - cost) * (1 - cost) * (1 - cost) * (1 - cost));
		};
		float rprob = SchlicksApproximation(r0, cost1);

		if (cost2 > 0.0f && RND2 > rprob)
			ray.dir = ((ray.dir * n) + (hitpointnrml *((n*cost1)*(n* cost1) - cost2))).Norm();
		else
			ray.dir = (ray.dir + hitpointnrml * (cost1 * 2.0f)).Norm();

		Vector3f tmpcolor;
		Trace(ray, scene, depth + 1, tmpcolor, hal, hal2);
		color = color + tmpcolor * 1.15f;
	}

	return;
}

int main() {
	Scene scene;
	Hitable* sphere;
	sphere = new Sphere(1.05f, Vector3f(-0.75f, -1.45f, -4.4f), Vector3f(6.0f, 6.0f, 6.f), 1.0f, Hitable::Type::Specular);
	scene.Add(sphere);
	scene.Add(new Sphere(0.6f, Vector3f(-1.75f, -1.95f, -3.1f), Vector3f(6.0f, 0.0f, 6.0f), 1.0f, Hitable::Type::Diffuse));
	scene.Add(new Sphere(0.7f, Vector3f(1.0f, -1.f, -3.7f), Vector3f(1.0f, 1.0f, 1.0f), 1.0f, Hitable::Type::Refraction));

	// "lights"
	scene.Add(new Sphere(0.5f, Vector3f(1.0f, 2.5f, -4.0f), Vector3f(10.0f, 10.0f, 10.0f), 1000, Hitable::Type::Diffuse));
	scene.Add(new Sphere(0.3f, Vector3f(-1.5f, 2.5f, -4.0f), Vector3f(10.f, 10.f, 10.0f), 10000, Hitable::Type::Diffuse));
		
	scene.Add(new Plane(2.5, Vector3f(0, 1, 0), Vector3f(6.f, 6.f, 6.f), 1, Hitable::Type::Diffuse)); // Bottom plane
	scene.Add(new Plane(2.5, Vector3f(0, -1, 0), Vector3f(6.f, 6.f, 6.f), 1, Hitable::Type::Diffuse)); // top plane
	scene.Add(new Plane(2.75, Vector3f(1, 0, 0), Vector3f(6.f, 0.f, 0.f), 1, Hitable::Type::Diffuse)); // left plane
	scene.Add(new Plane(2.5, Vector3f(-1, 0, 0), Vector3f(0.f, 6.f, 0.f), 1, Hitable::Type::Diffuse)); // right plane
	scene.Add(new Plane(5.5, Vector3f(0.0, 0, 1.0).Norm(), Vector3f(6.f, 6.f, 6.f), 1, Hitable::Type::Diffuse)); // back plane
	scene.Add(new Plane(-6.0, Vector3f(0.0, 0, -1.0).Norm(), Vector3f(0.f, 0.f, 0.f), 1, Hitable::Type::Diffuse)); // front plane

	std::vector<std::vector<Vector3f>> pix(Globals::width, std::vector<Vector3f>(Globals::width, Vector3f(0.f, 0.f, 0.f)));

	Halton hal, hal2;
	hal.Number(0, 2);
	hal2.Number(0, 2);

	const float spp = Globals::spp;

#pragma omp parallel for schedule(dynamic) firstprivate(hal,hal2)
	for (int col = 0; col < Globals::width; col++) {
		fprintf(stdout, "\rRendering: %1.0fspp %8.2f%%", spp, (float)col / Globals::width * 100.f);
		for (int row = 0; row < Globals::width; row++) {
			for (int s = 0; s < spp; s++) {
				Vector3f color;
				Ray ray;

				// rays start out from here
				ray.origin = Vector3f(0.f, 0.f, 0.f); 

				// construct image plane coordinates
				Vector3f cam = CamCords((float)col, (float) row); 

				// anti-aliasing
				cam.x = cam.x + RND / 700; 
				cam.y = cam.y + RND / 700;

				// point from the origin to the camera plane
				ray.dir = (cam - ray.origin).Norm(); 
				Trace(ray, scene, 0, color, hal, hal2);

				// write the contributions from every hit
				pix[col][row] = pix[col][row] + color / spp; 
			}
		}
	}

#pragma warning (disable : 4996)
	FILE *f = fopen(Globals::filename, "w");
	fprintf(f, "P3\n%d %d\n%d\n ", Globals::width, Globals::height, 255);
	for (int row = 0; row < Globals::height; row++) {
		for (int col = 0; col < Globals::width; col++) {
			fprintf(f, "%d %d %d ", std::min((int)pix[col][row].x, 255), std::min((int)pix[col][row].y, 255), std::min((int)pix[col][row].z, 255));
		}
		fprintf(f, "\n");
	}
	fclose(f);


	return 0;
}

