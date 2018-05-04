#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>		// MILO
#include "erand48.inc"	// MILO
#include <iostream>

#define M_PI 3.141592653589793238462643	// MILO

#define M_SAMPLES 4 // Karlex
#define USE_TENT_FILTER // Karlex

#define BYTE_RANGE 256

// targa file header 
typedef struct
{
	char id_length;
	char map_type;
	char image_type;
	int map_first;
	int map_length;
	char map_entry_size;
	int x;
	int y;
	int width;
	int height;
	char bits_per_pixel;
	char misc;
}targa_header;

inline int remainder(int number) { return number % BYTE_RANGE; }
inline int quotient(int number) { return number / BYTE_RANGE; }

struct Vec
{        // Usage: time ./smallpt 5000 && xv image.ppm
	double x, y, z;                  // position, also color (r,g,b)

	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec mult(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x*x + y * y + z * z)); }
	double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }

	// cross:
	Vec operator%(const Vec &b) { return Vec(y*b.z - z * b.y, z*b.x - x * b.z, x*b.y - y * b.x); }
};

struct Ray
{
	Vec o, d;
	Ray(const Vec &o_, const Vec &d_) : o(o_), d(d_) {}
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere
{
	double rad;       // radius
	Vec p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

	// returns distance, 0 if nohit
	double intersect(const Ray &r) const
	{ 
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) return 0; else det = sqrt(det);
		return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
	}
};

//Scene: radius, position, emission, color, material
Sphere spheres[] =
{
	Sphere(1e5, Vec(1e5 + 1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
	Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
	Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),Vec(),           DIFF),//Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
	Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
	Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
	Sphere(600, Vec(50,681.6 - .27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

// find hit index and distance if hit
inline bool intersect(const Ray &r, double &t, int &id)
{
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i--;)
		if ((d = spheres[i].intersect(r)) && d < t)
		{
			t = d;
			id = i;
		}
	return t < inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
	double t;                               // distance to intersection
	int id = 0;                               // id of intersected object
	if (!intersect(r, t, id)) return Vec(); // if miss, return black

	const Sphere &obj = spheres[id];        // the hit object
	Vec x = r.o + r.d*t;					// the hit pos
	Vec n = (x - obj.p).norm();				// outside normal
	Vec nl = n.dot(r.d) < 0 ? n : n * -1;	// radiance normal (outside normal or inside normal)
	Vec f = obj.c;							// sphere color
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl

	// Don¡¯t do Russian Roulette until after depth 5
	if (++depth > 5)
	{
		if (erand48(Xi) < p)
			f = f * (1 / p); // enlarge p to 1, lighten the color
		else return obj.e; //R.R.
	}

	if (depth > 100) return obj.e; // MILO

	if (obj.refl == DIFF)  // Ideal DIFFUSE reflection
	{                 
		double r1 = 2 * M_PI * erand48(Xi); // angle around
		double r2 = erand48(Xi); 
		double r2s = sqrt(r2); // distance from center
		Vec w = nl;
		Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u;
		Vec d = (u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2)).norm(); // {u, v, w} rotate matrix;
		return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
	}
	else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
		return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));

	Ray reflRay(x, r.d - n * 2 * n.dot(r.d));     // Ideal dielectric REFRACTION
	bool into = n.dot(nl) > 0;                // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
	if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0)    // Total internal reflection
		return obj.e + f.mult(radiance(reflRay, depth, Xi));
	Vec tdir = (r.d*nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0)*c*c*c*c*c, Tr = 1 - Re, P = .25 + .5*Re, RP = Re / P, TP = Tr / (1 - P);
	return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
		radiance(reflRay, depth, Xi)*RP : radiance(Ray(x, tdir), depth, Xi)*TP) :
		radiance(reflRay, depth, Xi)*Re + radiance(Ray(x, tdir), depth, Xi)*Tr);
}

// write targa file header
void write_header(targa_header header, FILE *tga)
{
	fputc(header.id_length, tga);
	fputc(header.map_type, tga);
	fputc(header.image_type, tga);

	fputc(remainder(header.map_first), tga);
	fputc(quotient(header.map_first), tga);

	fputc(remainder(header.map_length), tga);
	fputc(quotient(header.map_length), tga);

	fputc(header.map_entry_size, tga);

	fputc(remainder(header.x), tga);
	fputc(quotient(header.x), tga);
	fputc(remainder(header.y), tga);
	fputc(quotient(header.y), tga);

	fputc(remainder(header.width), tga);
	fputc(quotient(header.width), tga);
	fputc(remainder(header.height), tga);
	fputc(quotient(header.height), tga);

	fputc(header.bits_per_pixel, tga);
	fputc(header.misc, tga);
}

// write targa file
void write_tga(Vec *color, int width, int height)
{
	FILE *tga;
	targa_header header;

	int x, y;

	header.id_length = 0;
	header.map_type = 0;
	header.image_type = 2;

	header.map_first = 0;
	header.map_length = 0;
	header.map_entry_size = 0;

	header.x = 0;
	header.y = 0;
	header.width = width;
	header.height = height;

	header.bits_per_pixel = 24;
	header.misc = 0;

	tga = fopen("image.tga", "wb");
	write_header(header, tga);

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++)
		{
			// B G R order
			int i = (height - y - 1) * width + x;
			fputc(toInt(color[i].z), tga);
			fputc(toInt(color[i].y), tga);
			fputc(toInt(color[i].x), tga);
		}

	fclose(tga);
}

int main(int argc, char *argv[])
{
	clock_t start = clock(); // MILO
	int w = 256, h = 256; 
	int samps = argc == 2 ? atoi(argv[1]) / 4 : M_SAMPLES; // # samples, spp = samples * 4 (4 subpixel per pixel)
	Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
	Vec cx = Vec(w*.5135 / h);
	Vec cy = (cx%cam.d).norm()*.5135;
	Vec r;
	Vec *c = new Vec[w*h];

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
	for (int y = 0; y < h; y++)
	{                       // Loop over image rows
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.*y / (h - 1));
		unsigned short Xi[3] = { 0,0,y*y*y }; // MILO

		for (unsigned short x = 0; x < w; x++)   // Loop cols
		{
			for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++)     // 2x2 subpixel rows
			{
				for (int sx = 0; sx < 2; sx++, r = Vec())
				{        // 2x2 subpixel cols
					for (int s = 0; s < samps; s++)
					{
						// tent filter, nonlinear mapping r1 ¡Ê [0,2] to dx ¡Ê [-1,1]
						// to make sure sample points distribute denser on center£¬sparser on border
						// why do this? to make the simulation more accurate even few sample times?
						double r1 = 2 * erand48(Xi);
						double r2 = 2 * erand48(Xi);

						// seems no big difference
#ifdef USE_TENT_FILTER
						double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
#else
						double dx = r1 - 1;
						double dy = r2 - 1;
#endif

						Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
							cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;

						r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi)*(1. / samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior

					c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
				}
			}
		}
	}
	printf("\n%f sec\n", (float)(clock() - start) / CLOCKS_PER_SEC); // MILO

// 	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
// 
// 	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
// 
// 	for (int i = 0; i < w*h; i++)
// 		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));

	write_tga(c, w, h);

	system("pause"); // system("read")
}