#ifndef __SIMCONFIG_INCLUDE_GUARD__
#define __SIMCONFIG_INCLUDE_GUARD__
#include <string>

struct SimConfig
{
	int amount; //< amount of particles
	unsigned long long int simulation_steps; //< number of simulation steps
	
	unsigned width; //< x-coordinate, width of simulation area
	unsigned height; //< y-coordinate, width of simulation area
	unsigned depth; //< z-coordinate, width of simulation area
	
	double * x; //< arrays of x-coordinates of particles
	double * y; //< arrays of y-coordinates of particles
	double * z; //< arrays of z-coordinates of particles
	
	double * m; //< arrays of masses of particles
	
	double * vx; //< arrays of x-axis speeds of particles
	double * vy; //< arrays of y-axis speeds of particles 
	double * vz; //< arrays of z-axis speeds of particles
	
	int max_speed; //< maximum speed of particle
	int max_weight; //< maximum weight of particle
	
	std::string output_file; //< output_file, NOT USED
	
	SimConfig();
					
	~SimConfig();
	
};

#endif

