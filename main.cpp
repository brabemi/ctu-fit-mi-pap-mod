#include <iostream>

#define PARALLEL_OPENMP //< define for OpenMP, undefine for sequential
#include <omp.h> // OpenMP library

// if the following variable is NOT defined, program will not use any function nor include the CImg library
// this also means that other libraries are needed
#define CIMG_VISUAL

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
#define F_QUOC 0.0005 // compensatory quotient for  MAX_WEIGHT
#define BOUNCE_LOSS 0.8 // conversion rate of velocity on bounce, eg. 0.8 => 80% of speed after bounce

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
	bool graphics;

	// process input
	if(argc == 4)
	{
		threads = atoi(argv[1]);
		graphics = atoi(argv[2]) == 1;
		printf("Processing input file %s\n", argv[3]);
		processInputFile(argv[3], sconf);
	}
	else
	{
		printf("ERROR: Wrong/Missing parameters.\n");
		printf("Expect: \t%s THREADS GRAPHICS[0,1] INPUT_FILE", argv[0]);
		return 0;
	}
	
	printf("Input file processed.\n");	
	
	// iterators
	int i, j;
	// random gen. init
	srand(time(NULL));
	
	// pointers to arrays
	float * x = sconf.x;
	float * y = sconf.y;
	float * z = sconf.z;
	
	float * m = sconf.m;
	
	float * vx = sconf.vx;
	float * vy = sconf.vy;
	float * vz = sconf.vz;
	
	float * xnew = new float[sconf.amount];
	float * ynew = new float[sconf.amount];
	float * znew = new float[sconf.amount];
	
	// variables used in computations
	float ax, ay, az;
	float dx, dy, dz;
	float invr, invr3;
	float f;
	
	// definitino of "n" - used in algorithm on site
	int n = sconf.amount;
	
	// constants for computing particle movement 
	float dt = 0.1; // original value was 0.0001, that was too little for current values of particle parameters 
	float eps = 0.005;
	
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
			const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
			img.draw_point(x[i],y[i],color);
		}
		// create a Window (caption Playground) and fill it with image	
		main_disp = CImgDisplay(img,"Playground");
	}
	#endif

	unsigned steps = 0;

	float t1 = omp_get_wtime();

	printf("Starting simulation ...\n");
	while (steps < sconf.simulation_steps) {
		
		// compute new coordinates of all particles in parallel
		#ifdef PARALLEL_OPENMP
		#pragma omp parallel for num_threads(threads) private(ax,ay,az,dx,dy,dz,invr,invr3,f)
		#endif
		for(i=0; i<n; i++) { /* Foreach particle "i" ... */
			ax=0.0;
			ay=0.0;
			az=0.0;
			
			// reduction = parallel reduction of accelerations ax, ay, az
			#ifdef PARALLEL_OPENMP
			#pragma omp parallel for num_threads(threads) private(dx,dy,dz,invr,invr3,f) reduction(+:ax,ay,az)
			#endif
			for(j=0; j<n; j++) { /* Loop over all particles "j" */
				dx=x[j]-x[i];
				dy=y[j]-y[i];
				dz=z[j]-z[i];
				invr = 1.0/sqrtf(dx*dx + dy*dy + dz*dz + eps);
				//~ float tmp_sum = dx*dx + dy*dy + dz*dz + eps;
				//~ SSERsqrt(&invr, &tmp_sum);
				invr3 = invr*invr*invr;
				f=F_QUOC*m[j]*invr3;
				ax += f*dx; /* accumulate the acceleration from gravitational attraction */
				ay += f*dy;
				az += f*dz;
			}
			
			vx[i] += dt*ax; /* update velocity of particle "i" */
			vy[i] += dt*ay;
			vz[i] += dt*az;
			
			xnew[i] = x[i] + dt*vx[i] + 0.5*dt*dt*ax; /* update position of particle "i" */
			ynew[i] = y[i] + dt*vy[i] + 0.5*dt*dt*ay;
			znew[i] = z[i] + dt*vz[i] + 0.5*dt*dt*az;
			
			if( bounce(xnew[i], ynew[i], znew[i], sconf) ) {
				// update of particle velocity, change direction and value (BOUNCE_LOSS)
				vx[i] = debounce_vel(vx[i], xnew[i], 0, WIDTH);
				vy[i] = debounce_vel(vy[i], ynew[i], 0, HEIGHT);
				vz[i] = debounce_vel(vz[i], znew[i], 0, DEPTH);
				// update of particle position
				xnew[i] = debounce_pos(xnew[i], 0, WIDTH);
				ynew[i] = debounce_pos(ynew[i], 0, HEIGHT);
				znew[i] = debounce_pos(znew[i], 0, DEPTH);
			}
		}
		
		#ifdef CIMG_VISUAL
		if(graphics) {
			if(main_disp.is_closed()) break;
			
			for(i = 0; i < n; i++) {
				const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
				img.draw_circle(x[i],y[i],1,color);
				//~ img.draw_circle(x[i],y[i],2,green);
			}
		}
		#endif
		
		// update all coordinates of all particles in parallel
		#ifdef PARALLEL_OPENMP
		#pragma omp parallel for num_threads(threads)
		#endif
		for(i=0; i<n; i++) { /* copy updated positions back into original arrays */
			if( bounce(xnew[i], ynew[i], znew[i], sconf) ) {
				printf("Particle %d out of borders (x, y, z) = (%0.3f, %0.3f, %0.3f)\n", i, xnew[i], ynew[i], znew[i]);
			}
			x[i] = xnew[i];
			y[i] = ynew[i];
			z[i] = znew[i];
		}
			
		#ifdef CIMG_VISUAL
		if(graphics) {
			// redraw the image and show it in the window
			if(steps%500 == 0) {
				img.fill(0); //< black background 
			}
	
			// draw all particles
			for(i = 0; i < n; i++) {
				//~ const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
				img.draw_circle(x[i],y[i],1,red);
				//~ img.draw_circle(x[i],y[i],2,green);
			}
			// display the image in window	
			img.display(main_disp);
			
			// wait for some time
			cimg::wait(20); // in milisec
		}
		#endif
		if(100*steps % sconf.simulation_steps == 0) {
			printf("%.1f%% completed\n", 100.0*steps/sconf.simulation_steps);
		}
		steps++;
	}
	
	float t2 = omp_get_wtime(); // in seconds
	
	printf("Time: %f seconds\n",(t2-t1));
	
	delete [] xnew;
	delete [] ynew;
	delete [] znew;	

	return 0;
}
