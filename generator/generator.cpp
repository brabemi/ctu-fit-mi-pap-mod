#include <iostream>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* atoi */

#include <sstream> // stringstream

#include "SimConfig.h"

#include "ioproc.h"

using namespace std;

// Default configuration
void initDefaultSimConfig(SimConfig & sconf)
{
	sconf.width = 640;
	sconf.height = 640;
	sconf.depth = 500;
	
	sconf.amount = 100;
	sconf.simulation_steps = 1000;
	
	sconf.max_speed = 4;
	sconf.max_weight = 65536;
}

// Generates random particles within specified area in specified amount.
void generateRandomTestParticles(SimConfig & sconf)
{
	int i = 0;
	for(i = 0; i < sconf.amount; i++)
	{		
		sconf.x[i] = rand() % sconf.width;
		sconf.y[i] = rand() % sconf.height;
	    //~ sconf.z[i] = rand() % sconf.depth;
		sconf.z[i] = 0; // 2D simulation
		
		sconf.m[i] = rand() % sconf.max_weight;
		//~ sconf.vz[i] = rand() % 10;
		sconf.vx[i] = -1 * sconf.max_speed + rand() % (2 * sconf.max_speed + 1);
		sconf.vy[i] = -1 * sconf.max_speed + rand() % (2 * sconf.max_speed + 1);
		//~ sconf.vz[i] = -1 * sconf.max_speed + rand() % (2 * sconf.max_speed + 1);
		sconf.vz[i] = 0; // 2D simulation
	}
}

// Allocates dynamic arrays for particles.
int allocateDynamicArrays(SimConfig & sconf)
{
	if(sconf.amount < 1) return -1; // NOPE
	
	sconf.x = new double[sconf.amount];
	sconf.y = new double[sconf.amount];
	sconf.z = new double[sconf.amount];
	
	sconf.m = new double[sconf.amount];
	
	sconf.vx = new double[sconf.amount];
	sconf.vy = new double[sconf.amount];
	sconf.vz = new double[sconf.amount];
	
	return 0;
}

// Interactive command line setup of simulation.
void interactiveCommLineInput(SimConfig & sconf)
{
	cout << "Specify simulation parameters (Enter (only Enter!!!) to set the default value):" << endl;
	string inputString = "";
	
	// amount
	cout << endl << "Set amount of particles (default: " << sconf.amount << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.amount;
	} 
		
	// sim. steps
	cout << endl << "Set simulation steps (default: " << sconf.simulation_steps << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.simulation_steps;
	} 
	
	// width
	cout << endl << "Set sim. zone height (x-axis) (default: " << sconf.width << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.width;
	} 
	
	// height
	cout << endl << "Set sim. zone height (y-axis) (default: " << sconf.height << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.height;
	} 
	
	// depth
	cout << endl << "Set sim. zone depth (z-axis) (default: " << sconf.depth << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.depth;
	} 
	
	// max_speed
	cout << endl << "Set particle max. speed (default: " << sconf.max_speed << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.max_speed;
	} 
	
	// max_weight
	cout << endl << "Set particle max. weight (default: " << sconf.max_weight << "): ";
	getline(cin, inputString);
	if(inputString.length() != 0)
	{
		istringstream iss(inputString);
		iss >> sconf.max_weight;
	} 
		
}

// Main function - logic of generator
int main(int argc, char **argv)
{
	// create default configuration
	SimConfig sconf;
	
	// initialize default config
	initDefaultSimConfig(sconf);
	
	// default output filename
	const char * output_file = "output.txt";
	
	// process input - command line arguments
	if(argc == 1) // generates default config, random particles
	{
		// will use default values specified in function above
		cout << "No other parameters specified, will generate random particles with default setup and save it into file named \"" << output_file << "\"." << endl;
	}
	else if (argc == 9) // all arguments - run: "programname amount sim_steps width height depth max_speed max_weight output_file"
	{
		sconf.amount = atoi(argv[1]);
		sconf.simulation_steps = atoi(argv[2]); // possible failure - simulationSteps is unsigned long long int, atoi() returns only int
		
		sconf.width = atoi(argv[3]);
		sconf.height = atoi(argv[4]);
		sconf.depth = atoi(argv[5]);
		
		sconf.max_speed = atoi(argv[6]);
		sconf.max_weight = atoi(argv[7]);
		
		output_file = argv[8];
		
		cout << "Generating simulation settings and particles to file \"" << output_file << "\"..." << endl;
	}
	else if (argc == 2) // interactive setup of settings, run e.g.: "programname output_file"
	{
		interactiveCommLineInput(sconf);
		
		output_file = argv[1];
		
		cout << "Generating simulation settings and particles to file \"" << output_file << "\"..." << endl;
	}
	else
	{
		printf("GENERATOR ERROR: Weird amount of arguments, abort.");
	}
	
	if(allocateDynamicArrays(sconf)) // allocate arrays for particles
	{
		cout << "GENERATOR ERROR: Value of sconf.amount is less than 1." << endl;
		return 0;
	}
	
	generateRandomTestParticles(sconf); // fill the arrays
			
	createOutputFile(output_file,sconf); // write simulation details into output_file
	
	return 0;
}

