#include <iostream>
#include <omp.h> // OpenMP library

// if the following variable is NOT defined, program will not use any function nor include the CImg library
// this also means that other libraries are needed
#define CIMG_VISUAL

#ifdef CIMG_VISUAL
#include "CImg/CImg.h" // lib for visualisation
#else
#include <cstdio> // printf
#include <cstdlib> // srand
#include <cmath> // sqrt
#endif

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
- infinite simulation of NBody
- bouncing, borders
- no particle collisions
- no detection if 

TODO:
- define maximum steps of simulation
- input (definition of particles) from file/as parameters from comm. line
- input constants as parameters (WIDTH, MAX_SPEED, ...)
...
*/

#ifdef CIMG_VISUAL
using namespace cimg_library; // -> no need to use cimg_library::function()
#else
using namespace std;
#endif

bool bounce (double x, double y, double z, SimConfig & sconf) {
	return (x < 0) || (WIDTH < x) || (y < 0) || (HEIGHT < y) || (z < 0) || (DEPTH < z);
}

double debounce_vel (double vel, double pos, int min, int max) {
	if((pos < min) || (max < pos)) {
		return -1*BOUNCE_LOSS*vel;
	}
	return BOUNCE_LOSS*vel;
}

double debounce_pos (double pos, int min, int max) {
	if(pos < min) {
		return -1 * pos;
	}
	if(pos > max) {
		return max - (pos - max);
	}
	return pos;
}

int main(int argc, char** argv) {
	
	// simulation configuration
	SimConfig sconf;
	
	// process input
	if(argc == 2)
	{
		printf("Processing input file %s\n", argv[1]);
		processInputFile(argv[1], sconf);
	}
	else
	{
		printf("ERROR: Wrong/Missing parameters.\n");
		return 0;
	}
	
	printf("Input file processed.\n");
	
	#ifdef CIMG_VISUAL
	// image ~ "drawing panel"
	CImg<unsigned char> img( sconf.width, sconf.height,1,3);
	// colours
	const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
	#endif
	
	// iterators
	int i, j;
	// random gen. init
	srand(time(NULL));
	
	// pointers to arrays
	double * x = sconf.x;
	double * y = sconf.y;
	double * z = sconf.z;
	
	double * m = sconf.m;
	
	double * vx = sconf.vx;
	double * vy = sconf.vy;
	double * vz = sconf.vz;
	
	double * xnew = new double[sconf.amount];
	double * ynew = new double[sconf.amount];
	double * znew = new double[sconf.amount];
	
	// variables used in computations
	double ax, ay, az;
	double dx, dy, dz;
	double invr, invr3;
	double f;
	
	// definitino of "n" - used in algorithm on site
	int n = sconf.amount;
	
	// constants for computing particle movement 
	double dt = 0.1; // original value was 0.0001, that was too little for current values of particle parameters 
	double eps = 0.0001;
	
	#ifdef CIMG_VISUAL
	// initialization of window
	img.fill(0); //< fill img with black colour

	// draw all points
	for(i = 0; i < n; i++) {
		const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
		img.draw_point(x[i],y[i],color);
	}
	// create a Window (caption Playground) and fill it with image	
	CImgDisplay main_disp(img,"Playground");
	#endif

	unsigned steps = 0;
	
	double t1 = omp_get_wtime();
	printf("Starting simulation ...\n");
	while (steps < sconf.simulation_steps) {
		
		// compute new coordinates of all particles in parallel
		#pragma omp parallel for num_threads(4) private(ax,ay,az,dx,dy,dz,invr,invr3,f)
		for(i=0; i<n; i++) { /* Foreach particle "i" ... */
			ax=0.0;
			ay=0.0;
			az=0.0;
			
			// reduction = parallel reduction of accelerations ax, ay, az
			#pragma omp parallel for num_threads(4) private(dx,dy,dz,invr,invr3,f) reduction(+:ax,ay,az)
			for(j=0; j<n; j++) { /* Loop over all particles "j" */
				dx=x[j]-x[i];
				dy=y[j]-y[i];
				dz=z[j]-z[i];
				invr = 1.0/sqrt(dx*dx + dy*dy + dz*dz + eps);
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
		if(main_disp.is_closed()) break;
		
		for(i = 0; i < n; i++) {
			const unsigned char color[] = {(unsigned char) (m[i]/SQRT_MAX_WEIGHT), 255, (unsigned char) (((int) m[i])%SQRT_MAX_WEIGHT)};
			img.draw_circle(x[i],y[i],1,color);
			//~ img.draw_circle(x[i],y[i],2,green);
		}
		#endif
		
		// update all coordinates of all particles in parallel
		#pragma omp parallel for num_threads(4)
		for(i=0; i<n; i++) { /* copy updated positions back into original arrays */
			if( bounce(xnew[i], ynew[i], znew[i], sconf) ) {
				printf("Particle %d out of borders (x, y, z) = (%0.3f, %0.3f, %0.3f)\n", i, xnew[i], ynew[i], znew[i]);
			}
			x[i] = xnew[i];
			y[i] = ynew[i];
			z[i] = znew[i];
		}
			
		#ifdef CIMG_VISUAL
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
		#endif
		
		steps++;
	}
	
	double t2 = omp_get_wtime(); // in seconds
	
	printf("Time: %lf seconds\n",(t2-t1));
	
	delete [] xnew;
	delete [] ynew;
	delete [] znew;	

	return 0;
}
