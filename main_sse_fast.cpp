#include <iostream>

#define PARALLEL_OPENMP //< define for OpenMP, undefine for sequential
#include <omp.h> // OpenMP library

// if the following variable is NOT defined, program will not use any function nor include the CImg library
// this also means that other libraries are needed
//~ #define CIMG_VISUAL

#define SSE_SQRT

#define LOGGING // all logging output except for the line with "#THREADS #SECONDS" info
#undef LOGGING

#ifdef CIMG_VISUAL
#include "CImg/CImg.h" // lib for visualisation
#endif
#include <cstdio> // printf
#include <cstdlib> // srand
#include <cmath> // sqrt
#include <xmmintrin.h> //SSERsqrt

// include our functions/structs
#include "generator/ioproc.h" // process input file
#include "generator/SimConfig.h" // Struct with simulation config/settings

// current mapping of simulation parameters - just to clarify what is what 
#define WIDTH sconf.width // // x-coordinate, width of simulation area
#define HEIGHT sconf.height // y-coordinate, height of simulation area
#define DEPTH sconf.depth // z-coordinate, depth of simulation area
#define MAX_SPEED sconf.max_speed // maximal default speed of particle, eg. 4 => speed in (-4;4)
#define MAX_WEIGHT sconf.max_weight // maximal mass of particle, eg. 4 => speed in (0;4)
// amount of particles
#define AMOUNT sconf.amount // amount of particles
// number of simulation steps is specified in SimConfig struct as simulation_steps
#define SQRT_MAX_WEIGHT (int)sqrt(sconf.max_weight) // constant for calculation particle color {mass/SQRT_MAX_WEIGHT, 255, mass%SQRT_MAX_WEIGHT}
// ------------- END OF MODIFIED MAPPINGs

// parameters for calculations (movement of particles)
#define F_QUOC 0.0005f // compensatory quotient for  MAX_WEIGHT
#define BOUNCE_LOSS 0.8f // conversion rate of velocity on bounce, eg. 0.8 => 80% of speed after bounce

/*
CURRENT STATE:
- finite simulation of NBody
- bouncing, borders
- no particle collisions
- no detection if 

LAST UPDATE (just to clarify what happened - move to CURRENT STATE after review):
- inner loop in simulation is now parallelized with parallel reduction of accelerations
- define maximum steps of simulation - viz SimConfig::simulation_steps
- input (definition of particles) from file/as parameters from comm. line - viz generator.cpp and ioproc.cpp
- input constants as parameters (WIDTH, MAX_SPEED, ...) - viz generator.cpp
- for CImg visualisation, variable "CIMG_VISUAL" must be defined - uncomment #define CIMG_VISUAL at the beginning of this file

- created standalone generator (in directory "generator")
- simulation should be run with: ./simulator input.txt
- valgrind says the program leaks a bit (6 block every time) - with --leak-check=full --show-leak-kinds=all, it seems that libgomp causes that

- reworked makefile (make compile, make clean)

TODO:
...

*/

#ifdef CIMG_VISUAL
using namespace cimg_library; // -> no need to use cimg_library::function()
#else
using namespace std;
#endif

typedef union
{
  __m128 * v;
  float * f;
} f4chunk;

typedef union
{
  __m128 v;
  float f [4];
} f4vector;

bool bounce (float x, float y, float z, SimConfig & sconf) {
	return (x < 0) || (WIDTH < x) || (y < 0) || (HEIGHT < y) || (z < 0) || (DEPTH < z);
}

float debounce_vel (float vel, float pos, int min, int max) {
	if((pos < min) || (max < pos)) {
		return -1*BOUNCE_LOSS*vel;
	}
	return BOUNCE_LOSS*vel;
}

float debounce_pos (float pos, int min, int max) {
	if(pos < min) {
		return -1 * pos;
	}
	if(pos > max) {
		return max - (pos - max);
	}
	return pos;
}

inline void SSERsqrt( float * pOut, float * pIn )
{
	_mm_store_ss( pOut, _mm_rsqrt_ss( _mm_load_ss( pIn ) ) );
	// compiles to movss, sqrtss, movss
}


int main(int argc, char** argv) {
	
	// simulation configuration
	SimConfig sconf;
	int threads;

	#ifdef CIMG_VISUAL
	bool graphics;
	#endif

	// process input
	if(argc == 4)
	{
		threads = atoi(argv[1]);
		#ifdef CIMG_VISUAL
		graphics = atoi(argv[2]) == 1;
		#endif
		#ifdef LOGGING
		printf("Processing input file %s\n", argv[3]);
		#endif
		processInputFile(argv[3], sconf);
	}
	else
	{
		printf("ERROR: Wrong/Missing parameters.\n");
		printf("Expect: \t%s THREADS GRAPHICS[0,1] INPUT_FILE", argv[0]);
		return 0;
	}
	#ifdef LOGGING
	printf("Input file processed.\n");	
	#endif
	// iterators
	int i, j;
	// random gen. init
	srand(time(NULL));
	
	// pointers to arrays
	f4chunk x, y, z, m, vx, vy, vz, xnew, ynew, znew;
	x.f = sconf.x;
	y.f = sconf.y;
	z.f = sconf.z;
	//~ float * x = sconf.x;
	//~ float * y = sconf.y;
	//~ float * z = sconf.z;
	
	m.f = sconf.m;
	//~ float * m = sconf.m;
	
	vx.f = sconf.vx;
	vy.f = sconf.vy;
	vz.f = sconf.vz;
	//~ float * vx = sconf.vx;
	//~ float * vy = sconf.vy;
	//~ float * vz = sconf.vz;
	
	xnew.f = new float[sconf.amount];
	ynew.f = new float[sconf.amount];
	znew.f = new float[sconf.amount];
	//~ float * xnew = new float[sconf.amount];
	//~ float * ynew = new float[sconf.amount];
	//~ float * znew = new float[sconf.amount];
	
	// variables used in computations
	f4vector ax, ay, az;
	f4vector dx, dy, dz;
	f4vector invr, invr3;
	f4vector f;
	f4vector f_q;
	f_q.f[0] = f_q.f[1] = f_q.f[2] = f_q.f[3] = F_QUOC;
	//~ float ax, ay, az;
	//~ float dx, dy, dz;
	//~ float invr, invr3;
	//~ float f;
	
	// definitino of "n" - used in algorithm on site
	int n = sconf.amount;
	
	// constants for computing particle movement 
	f4vector dt; // original value was 0.0001, that was too little for current values of particle parameters 
	dt.f[0] = 0.1f; // original value was 0.0001, that was too little for current values of particle parameters 
	dt.f[1] = 0.1f; // original value was 0.0001, that was too little for current values of particle parameters 
	dt.f[2] = 0.1f; // original value was 0.0001, that was too little for current values of particle parameters 
	dt.f[3] = 0.1f; // original value was 0.0001, that was too little for current values of particle parameters 
	//~ float dt = 0.1; // original value was 0.0001, that was too little for current values of particle parameters 
	f4vector eps;
	eps.f[0] = 0.005f;
	eps.f[1] = 0.005f;
	eps.f[2] = 0.005f;
	eps.f[3] = 0.005f;
	//~ float eps = 0.005;
	
	#ifdef CIMG_VISUAL
	// image ~ "drawing panel"
	CImg<unsigned char> img;
	CImgDisplay main_disp;
	// colours
	const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };

	if(graphics) {
		// image ~ "drawing panel"
		img = CImg<unsigned char> ( sconf.width, sconf.height,1,3);

		// initialization of window
		img.fill(0); //< fill img with black colour
	
		// draw all points
		for(i = 0; i < n; i++) {
			const unsigned char color[] = {(unsigned char) (m.f[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m.f[i])%SQRT_MAX_WEIGHT)};
			img.draw_point(x.f[i],y.f[i],color);
		}
		// create a Window (caption Playground) and fill it with image	
		main_disp = CImgDisplay(img,"Playground");
	}
	#endif

	unsigned steps = 0;

	float t1 = omp_get_wtime();
	#ifdef LOGGING
	printf("Starting simulation ...\n");
	#endif
	while (steps < sconf.simulation_steps) {
		
		// compute new coordinates of all particles in parallel
		#ifdef PARALLEL_OPENMP
		#pragma omp parallel for num_threads(threads) private(ax,ay,az,dx,dy,dz,invr,invr3,f,j)
		#endif
		for(i=0; i<n; i++) { /* Foreach particle "i" ... */
			ax.f[0]=0.0f;
			ax.f[1]=0.0f;
			ax.f[2]=0.0f;
			ax.f[3]=0.0f;
			ay.f[0]=0.0f;
			ay.f[1]=0.0f;
			ay.f[2]=0.0f;
			ay.f[3]=0.0f;
			az.f[0]=0.0f;
			az.f[1]=0.0f;
			az.f[2]=0.0f;
			az.f[3]=0.0f;

			f4vector tmp_x, tmp_y, tmp_z;
			tmp_x.f[0] = tmp_x.f[1] = tmp_x.f[2] = tmp_x.f[3] = x.f[i];
			tmp_y.f[0] = tmp_y.f[1] = tmp_y.f[2] = tmp_y.f[3] = y.f[i];
			tmp_z.f[0] = tmp_z.f[1] = tmp_z.f[2] = tmp_z.f[3] = z.f[i];
			// reduction = parallel reduction of accelerations ax, ay, az
			//~ #ifdef PARALLEL_OPENMP
			//~ #pragma omp parallel for num_threads(threads) private(dx,dy,dz,invr,invr3,f) reduction(+:ax,ay,az)
			//~ #endif
			for(j=0; j<n/4; j++) { /* Loop over all particles "j" */
				dx.v=x.v[j]-tmp_x.v;
				dy.v=y.v[j]-tmp_y.v;
				dz.v=z.v[j]-tmp_z.v;

				//~ dx=x[j]-x[i];
				//~ dy=y[j]-y[i];
				//~ dz=z[j]-z[i];
				
				invr.v = _mm_rsqrt_ps(dx.v*dx.v + dy.v*dy.v + dz.v*dz.v + eps.v);
				
				
				//~ #ifdef SSE_SQRT
				//~ float tmp_sum = dx*dx + dy*dy + dz*dz + eps;
				//~ SSERsqrt(&invr, &tmp_sum);
				//~ #else
				//~ invr = 1.0/sqrtf(dx*dx + dy*dy + dz*dz + eps);
				//~ #endif
				
				invr3.v = invr.v*invr.v*invr.v;
				//~ invr3 = invr*invr*invr;
				f.v = f_q.v * m.v[j] * invr3.v;
				ax.v += f.v*dx.v; /* accumulate the acceleration from gravitational attraction */
				ay.v += f.v*dy.v;
				az.v += f.v*dz.v;
			}
			
			float tmp_ax = (ax.f[0] + ax.f[1] + ax.f[2] + ax.f[3]);
			float tmp_ay = (ay.f[0] + ay.f[1] + ay.f[2] + ay.f[3]);
			float tmp_az = (az.f[0] + az.f[1] + az.f[2] + az.f[3]);
			vx.f[i] += dt.f[0]*tmp_ax; /* update velocity of particle "i" */
			vy.f[i] += dt.f[0]*tmp_ay; /* update velocity of particle "i" */
			vz.f[i] += dt.f[0]*tmp_az; /* update velocity of particle "i" */
			
			xnew.f[i] = x.f[i] + dt.f[0]*vx.f[i] + 0.5f*dt.f[0]*dt.f[0]*tmp_ax; /* update position of particle "i" */
			ynew.f[i] = y.f[i] + dt.f[0]*vy.f[i] + 0.5f*dt.f[0]*dt.f[0]*tmp_ay; /* update position of particle "i" */
			znew.f[i] = z.f[i] + dt.f[0]*vz.f[i] + 0.5f*dt.f[0]*dt.f[0]*tmp_az; /* update position of particle "i" */
			
			if( bounce(xnew.f[i], ynew.f[i], znew.f[i], sconf) ) {
				// update of particle velocity, change direction and value (BOUNCE_LOSS)
				vx.f[i] = debounce_vel(vx.f[i], xnew.f[i], 0, WIDTH);
				vy.f[i] = debounce_vel(vy.f[i], ynew.f[i], 0, HEIGHT);
				vz.f[i] = debounce_vel(vz.f[i], znew.f[i], 0, DEPTH);
				// update of particle position
				xnew.f[i] = debounce_pos(xnew.f[i], 0, WIDTH);
				ynew.f[i] = debounce_pos(ynew.f[i], 0, HEIGHT);
				znew.f[i] = debounce_pos(znew.f[i], 0, DEPTH);
			}
		}
		
		#ifdef CIMG_VISUAL
		if(graphics) {
			if(main_disp.is_closed()) break;
			
			for(i = 0; i < n; i++) {
				const unsigned char color[] = {(unsigned char) (m.f[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m.f[i])%SQRT_MAX_WEIGHT)};
				img.draw_circle(x.f[i],y.f[i],1,color);
				//~ img.draw_circle(x.f[i],y.f[i],2,green);
			}
		}
		#endif
		
		// update all coordinates of all particles in parallel
		#ifdef PARALLEL_OPENMP
		#pragma omp parallel for num_threads(threads)
		#endif
		for(i=0; i<n; i++) { /* copy updated positions back into original arrays */
			if( bounce(xnew.f[i], ynew.f[i], znew.f[i], sconf) ) {
				#ifdef LOGGING
				printf("Particle %d out of borders (x, y, z) = (%0.3f, %0.3f, %0.3f)\n", i, xnew.f[i], ynew.f[i], znew.f[i]);
				#endif
			}
			x.f[i] = xnew.f[i];
			y.f[i] = ynew.f[i];
			z.f[i] = znew.f[i];
		}
			
		#ifdef CIMG_VISUAL
		if(graphics) {
			// redraw the image and show it in the window
			if(steps%500 == 0) {
				img.fill(0); //< black background 
			}
	
			// draw all particles
			for(i = 0; i < n; i++) {
				//~ const unsigned char color[] = {(unsigned char) (m.f[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m.f[i])%SQRT_MAX_WEIGHT)};
				img.draw_circle(x.f[i],y.f[i],1,red);
				//~ img.draw_circle(x.f[i],y.f[i],2,green);
			}
			// display the image in window	
			img.display(main_disp);
			
			// wait for some time
			cimg::wait(20); // in milisec
		}
		#endif
		if(100*steps % sconf.simulation_steps == 0) {
			#ifdef LOGGING
			printf("%.1f%% completed\n", 100.0*steps/((int)sconf.simulation_steps));
			#endif
		}
		steps++;
	}
	
	float t2 = omp_get_wtime(); // in seconds
	#ifdef LOGGING
	printf("Time: %f seconds\n",(t2-t1));
	#endif
	printf("%d %f\n", threads, (t2-t1));
	
	delete [] xnew.f;
	delete [] ynew.f;
	delete [] znew.f;	

	return 0;
}
