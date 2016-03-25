#include <iostream>

#include "SimConfig.h"

SimConfig::SimConfig(): amount(100), simulation_steps(1000), width(640), height(640), depth(500),
				x(NULL), y(NULL), z(NULL), m(NULL), vx(NULL), vy(NULL), vz(NULL),
				max_speed(4), max_weight(65536), output_file("output.txt") {}
					
SimConfig::~SimConfig()
{
	if(x != NULL) delete [] x;
	if(y != NULL) delete [] y;
	if(z != NULL) delete [] z;
		
	if(m != NULL) delete [] m;
		
	if(vx != NULL) delete [] vx;
	if(vy != NULL) delete [] vy;
	if(vz != NULL) delete [] vz;
}
