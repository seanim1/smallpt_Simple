// Simplified smallpt: https://raytracey.blogspot.com/2015/10/gpu-path-tracing-tutorial-1-drawing.html
// C++ smallpt: https://github.com/codesavory/very-own-path-tracer
// Kevin Beason smallpt Slide: https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Fwww.kevinbeason.com%2Fsmallpt%2Fsmallpt_cline_slides.ppt&wdOrigin=BROWSELINK
// based on smallpt, a Path Tracer by Kevin Beason, 2008

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>
#include "glm/glm.hpp"

#define M_PI 3.14159265359f  // pi
#define width 512  // screenwidth
#define height 384 // screenheight
#define samps 4 // samples 
#define TEMPORAL 0.75f
//typedef glm::vec3 vec3;
using namespace glm;

float getrandom()
{
	return (float)rand() / (float)RAND_MAX;
}

// Ray Structure
struct Ray
{
	vec3 orig; // ray origin
	vec3 dir;  // ray direction

	Ray(vec3 o_, vec3 d_) : orig(o_), dir(d_) {}
};

// Enum of material types used in RADIANCE function
enum Refl_t { DIFF, SPEC, REFR };

//smallpt only supports spheres, Sphere Structure
struct Sphere
{
	float rad;            // radius 
	vec3 pos, emi, col; // position, emission, colour 
	Refl_t refl;          // reflection type (e.g. diffuse)

	//returns distance or 0 if no hit
	float intersect_sphere(const Ray& r) const
	{
		vec3 op = pos - r.orig;    // distance from ray.orig to center sphere 
		float t, epsilon = 0.0001f;  // epsilon required to prevent floating point precision artefacts
		float b = dot(op, r.dir);    // b in quadratic equation
		float disc = b * b - dot(op, op) + rad * rad;  // discriminant quadratic equation
		if (disc < 0) return 0;       // if disc < 0, no real solution (we're not interested in complex roots) 
		else disc = sqrtf(disc);    // if disc >= 0, check for solutions using negative and positive discriminant
		return (t = b - disc) > epsilon ? t : ((t = b + disc) > epsilon ? t : 0); // pick closest point in front of ray origin
	}
};

// HARD CODED scene description
// The Scene Description consists of a bunch of spheres
// Scene: radius, position, emission, color, material
Sphere spheres[] = {
	{ 1e5f,		{ 1e5f + 1.0f, 40.8f, 81.6f },		{ 0.0f, 0.0f, 0.0f }, { 0.75f, 0.25f, 0.25f },	DIFF }, //Left 
	{ 1e5f,		{ -1e5f + 99.0f, 40.8f, 81.6f },	{ 0.0f, 0.0f, 0.0f }, { .25f, .25f, .75f },		DIFF }, //Rght 
	{ 1e5f,		{ 50.0f, 40.8f, 1e5f },				{ 0.0f, 0.0f, 0.0f }, { .75f, .75f, .75f },		DIFF }, //Back 
	{ 1e5f,		{ 50.0f, 40.8f, -1e5f + 600.0f },	{ 0.0f, 0.0f, 0.0f }, { 1.00f, 1.00f, 1.00f },	DIFF }, //Frnt 
	{ 1e5f,		{ 50.0f, 1e5f, 81.6f },				{ 0.0f, 0.0f, 0.0f }, { .75f, .75f, .75f },		DIFF }, //Botm 
	{ 1e5f,		{ 50.0f, -1e5f + 81.6f, 81.6f },	{ 0.0f, 0.0f, 0.0f }, { .75f, .75f, .75f },		DIFF }, //Top 
	{ 16.5f,	{ 27.0f, 16.5f, 47.0f },			{ 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f },		DIFF }, // small sphere 1
	{ 16.5f,	{ 73.0f, 16.5f, 78.0f },			{ 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f },		DIFF }, // small sphere 2
	{ 600.0f,	{ 50.0f, 681.6f - .77f, 81.6f },	{ 2.0f, 1.8f, 1.6f }, { 0.0f, 0.0f, 0.0f },		DIFF }  // Light
};

// CLAMP FUNCTION
inline float clamp(float x) { return x < 0.0f ? 0.0f : x > 1.0f ? 1.0f : x; }

// CONVERTS FLOATS TO INTEGERS TO BE SAVED in PPM File
inline int toInt(float x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }  // convert RGB float in range [0,1] to int in range [0, 255] and perform gamma correction of 2.2

// INTERSECTS ray with SCENE
inline bool intersect_scene(const Ray& r, float& t, int& id) {
	float n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;  // t is distance to closest intersection, initialise t to a huge number outside scene
	for (int i = int(n); i--;)  // test all scene objects for intersection
		if ((d = spheres[i].intersect_sphere(r)) && d < t) {  // if newly computed intersection distance d is smaller than current closest intersection distance
			t = d;  // keep track of distance along ray to closest intersection point 
			id = i; // and closest intersected object
		}
	return t < inf; // returns true if an intersection with the scene occurred, false when no hit
}

// COMPUTES the radiance estimate along the Ray r
vec3 radiance(Ray& r)
{
	vec3 accucolor = vec3(0.0f, 0.0f, 0.0f); // accumulates ray colour with each iteration through bounce loop
	vec3 mask = vec3(1.0f, 1.0f, 1.0f);

	// ray bounce loop (no Russian Roulette used) 
	for (int bounces = 0; bounces < 4; bounces++) {  // iteration up to 4 bounces (replaces recursion in CPU code)

		float t;           // distance to closest intersection 
		int id = 0;        // index of closest intersected sphere 

		// test ray for intersection with scene
		if (!intersect_scene(r, t, id))
			return vec3(0.0f, 0.0f, 0.0f); // if miss, return black

		// else, we've got a hit!
		// compute hitpoint and normal
		const Sphere& obj = spheres[id];  // hitobject
		vec3 x = r.orig + r.dir * t;          // hitpoint 
		vec3 n = normalize(x - obj.pos);    // normal
		vec3 nl = dot(n, r.dir) < 0 ? n : -n; // front facing normal

		// add emission of current sphere to accumulated colour
		// (first term in rendering equation sum) 
		accucolor += mask * obj.emi;

		// all spheres in the scene are diffuse
		// diffuse material reflects light uniformly in all directions
		// generate new diffuse ray:
		// origin = hitpoint of previous ray in path
		// random direction in hemisphere above hitpoint (see "Realistic Ray Tracing", P. Shirley)

		// create 2 random numbers
		float r1 = 2 * M_PI * getrandom(); // pick random number on unit circle (radius = 1, circumference = 2*Pi) for azimuth
		float r2 = getrandom();  // pick random number for elevation
		float r2s = sqrtf(r2);

		// compute local orthonormal basis uvw at hitpoint to use for calculation random ray direction 
		// first vector = normal at hitpoint, second vector is orthogonal to first, third vector is orthogonal to first two vectors
		vec3 w = nl;
		vec3 u = normalize(cross((fabs(w.x) > .1 ? vec3(0, 1, 0) : vec3(1, 0, 0)), w));
		vec3 v = cross(w, u);

		// compute random ray direction on hemisphere using polar coordinates
		// cosine weighted importance sampling (favours ray directions closer to normal direction)
		vec3 d = normalize(u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrtf(1 - r2));

		// new ray origin is intersection point of previous ray with scene
		r.orig = x + nl * 0.05f; // offset ray origin slightly to prevent self intersection
		r.dir = d;

		mask *= obj.col;    // multiply with colour of object       
		mask *= dot(d, nl);  // weigh light contribution using cosine of angle between incident light and normal
		mask *= 2;          // fudge factor
	}
	return accucolor;
}

int main(int argc, char* argv[])
{
	vec3* output = new vec3[width * height];										// The image matrix
	Ray cam(vec3(50, 52, 295.6), normalize(vec3(0.0, -0.042612, -1))); // first hardcoded camera ray(origin, direction)
	vec3 cx = vec3(width * .5135 / height, 0.0f, 0.0f); // ray direction offset in x direction
	vec3 cy = normalize(cross(cx, cam.dir)) * .5135f; // ray direction offset in y direction (.5135 is field of view angle)
	vec3 r = vec3(0); // r is final pixel color 

	auto start = std::chrono::high_resolution_clock::now();
	//Loop over all image pixels
	for (int y = 0; y < height; y++)	//loop over image rows
	{
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps, 100. * y / (height - 1));	// print progress
		for (unsigned short x = 0; x < width; x++)	// Loop Columns
		{
			unsigned int i = (height - y - 1) * width + x; // index of current pixel
			r = vec3(0.f);
			for (int s = 0; s < samps; s++)	// Camera rays are pushed ^^^^^ forward to start in interior
			{
				vec3 d = cam.dir + cx * ((.2f + x) / width - .5f) + cy * ((.2f + y) / height - .5f);
				Ray ray = Ray(cam.orig + d * 40.f, normalize(d));
				vec3 result = radiance(ray) * (1.f / samps);
				r = r + result;
			}
			// write rgb value of pixel to image buffer on the GPU, clamp value to [0.0f, 1.0f] range
			output[i] = vec3(clamp(r.x, 0.0f, 1.0f), clamp(r.y, 0.0f, 1.0f), clamp(r.z, 0.0f, 1.0f));
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> elapsed = end - start;
	std::cout << "\nElapsed time: " << elapsed.count() << " seconds." << std::endl;

	// Write out the file to a PPM
	FILE* f = fopen("smallpt_Simple.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
	vec3 gamma_correction = vec3(1 / 2.2);
	for (int i = 0; i < width * height; i++) {
		ivec3 pixel = pow(clamp(output[i], 0.f, 1.f), gamma_correction) * 255.f + .5f;
		fprintf(f, "%d %d %d ", pixel.x, pixel.y, pixel.z);
	}
}