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
	
	float * x; //< arrays of x-coordinates of particles
	float * y; //< arrays of y-coordinates of particles
	float * z; //< arrays of z-coordinates of particles
	
	float * m; //< arrays of masses of particles
	
	float * vx; //< arrays of x-axis speeds of particles
	float * vy; //< arrays of y-axis speeds of particles 
	float * vz; //< arrays of z-axis speeds of particles
	
	int max_speed; //< maximum speed of particle
	int max_weight; //< maximum weight of particle
	
	std::string output_file; //< output_file, NOT USED
	
	SimConfig();
					
	~SimConfig();
	
};

#endif

