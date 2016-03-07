#include <iostream>

#include <omp.h> // OpenMP library

#include "CImg/CImg.h" // lib for visualisation

// definition of parameters
#define WIDTH 1200 // // x-coordinate, width of simulation area
#define DEPTH 640 // y-coordinate, depth of simulation area
#define HEIGHT 640 // z-coordinate, height of simulation area
#define AMOUNT 500 // amount of particles

/*
CURRENT STATE:
- infinite simulation of NBody
- no bouncing, no borders
...

TODO:
- define maximum steps of simulation
- input (definition of particles) from file/as parameters from comm. line
...
*/

using namespace cimg_library; // -> no need to use cimg_library::function()

int main(int argc, char** argv) {
	
	// image ~ "drawing panel"
	CImg<unsigned char> img( WIDTH , HEIGHT ,1,3);
	// colours
	const unsigned char red[] = { 255,0,0 }, green[] = { 0,255,0 }, blue[] = { 0,0,255 };
	
	// iterators
	int i, j;
	// random gen. init
	srand(time(NULL));
	
	// data structures for particles
	
	// XYZ-coordinates
	int x[AMOUNT];
	int y[AMOUNT];
	int z[AMOUNT];

	int xnew[AMOUNT];
	int ynew[AMOUNT];
	int znew[AMOUNT];
	
	// masses
	int m[AMOUNT];
	
	// speeds
	int vx[AMOUNT];
	int vy[AMOUNT];
	int vz[AMOUNT];
	
	// generation of random particle parameter values
	for(i = 0; i < AMOUNT; i++)
	{
		x[i] = rand() % WIDTH;
		y[i] = rand() % HEIGHT;
		z[i] = rand() % DEPTH;
		
		m[i] = rand() % 100;
		
		vx[i] = rand() % 10;
		vy[i] = (rand() % 10);
		vz[i] = rand() % 10;
	}
	

	/*
		int x[] = {10,20,30,40,50,60};
		int y[] = {10,20,30,40,50,60};
		int z[] = {10,20,30,40,50,60};

		int xnew[] = {1,2,3,4,5,6};
		int ynew[] = {1,2,3,4,5,6};
		int znew[] = {1,2,3,4,5,6};

		int m[] = {10,20,30,40,50,60};

		int vx[] = {10,10,10,10,10,10};
		int vy[] = {1,2,3,4,5,6};
		int vz[] = {1,2,3,4,5,6};
	*/
	
	// variables used in computations
	double ax, ay, az;
	int dx, dy, dz;	
	double invr, invr3;
	double f;
	
	// definitino of "n" - used in algorithm on site
	int n = AMOUNT;
	
	// constants for computing particle movement 
	float dt = 0.1; // original value was 0.0001, that was too little for current values of particle parameters 
	float eps = 0.0001;
	
	// initialization of window
	img.fill(0); //< fill img with black colour

	// draw all points
	for(i = 0; i < n; i++)
		img.draw_point(x[i],y[i],green);
		
	// create a Window (caption Playground) and fill it with image	
	CImgDisplay main_disp(img,"Playground");

	// variables for detection of all particles in one place (point)
	int xFin, yFin, zFin;

	// while the Window is still opened ...
	while (!main_disp.is_closed()) {
		
		// compute new coordinates of all particles in parallel
		#pragma omp parallel for num_threads(4) private(ax,ay,az,dx,dy,dz,invr,invr3,f)
		for(i=0; i<n; i++) { /* Foreach particle "i" ... */
			ax=0.0;
			ay=0.0;
			az=0.0;

			for(j=0; j<n; j++) { /* Loop over all particles "j" */
				dx=x[j]-x[i];
				dy=y[j]-y[i];
				dz=z[j]-z[i];
				invr = 1.0/sqrt(dx*dx + dy*dy + dz*dz + eps);
				invr3 = invr*invr*invr;
				f=m[j]*invr3;
				ax += f*dx; /* accumulate the acceleration from gravitational attraction */
				ay += f*dy;
				az += f*dx;
			}

			xnew[i] = x[i] + dt*vx[i] + 0.5*dt*dt*ax; /* update position of particle "i" */
			ynew[i] = y[i] + dt*vy[i] + 0.5*dt*dt*ay;
			znew[i] = z[i] + dt*vz[i] + 0.5*dt*dt*az;
			vx[i] += dt*ax; /* update velocity of particle "i" */
			vy[i] += dt*ay;
			vz[i] += dt*az;
		}
		
		// update all coordinates of all particles in parallel
		#pragma omp parallel for num_threads(4)
		for(i=0; i<n; i++) { /* copy updated positions back into original arrays */
			x[i] = xnew[i];
			y[i] = ynew[i];
			z[i] = znew[i];
		}
		
		
		// indicate whether all particles are in the same place (same XYZ coordinates)
		xFin = x[0];
		yFin = y[0];
		zFin = z[0];
		
		bool same_place = true;
		
		for(i = 1; i < n; i++)
		{
			if((xFin == x[i]) &&
			(yFin == y[i])
			)
			{}
			else same_place = false;
		}
		
		// redraw the image and show it in the window
		img.fill(0); //< black background 
		
		// information about same coordinates of all particles
		if(same_place)
			img.draw_text(WIDTH/10, HEIGHT/2, "ALL IN ONE PLACE", red, blue);
			
		// draw all particles
		for(i = 0; i < n; i++)
			img.draw_circle(x[i],y[i],2,green);
			
		// display the image in window	
		img.display(main_disp);
		
		// wait for some time
		cimg::wait(20); // in milisec
	}

	return 0;
}