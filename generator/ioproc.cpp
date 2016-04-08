#include "ioproc.h"
#include "SimConfig.h"

#include <iostream> // cout, etc.
#include <fstream> // file
#include <sstream> // stringstream

using namespace std;

//#define DEBUG_INPUT // debug output (to cout) of input process function
//#define DEBUG_TEST_IO // run test (compile with main)

/**
 * Process input file - read particles from file.
 */
int processInputFile(const char * inputFilename, SimConfig & sconf)
{	
	// open file
	ifstream inputFile(inputFilename, ios::in);
	
	// check if the file is opened
	if(!inputFile.is_open())
	{
		cout << "ERROR: Could not open the input file." << endl;
		return 0;
	}
	
	// expected format of file:
	/* (first line of file): amount sim_steps width height depth max_speed max_weight
	 * (every other line): x y z m vx vy vz<lf> (x y (...) vz = floats)
	 */
	
	string inputLine = "";
	
	// read the header - amount of particles in file
	getline(inputFile,inputLine);
	
	// create stringstream and read the amount (?and other possible parameters?)
	istringstream issHeader(inputLine);
	issHeader >> sconf.amount;
	issHeader >> sconf.simulation_steps;
	
	issHeader >> sconf.width;
	issHeader >> sconf.height;
	issHeader >> sconf.depth;
	
	issHeader >> sconf.max_speed;
	issHeader >> sconf.max_weight;
	
	
	// check value of amount - should be more than 0
	if(sconf.amount < 1) 
	{
		cout << "ERROR: Amount of particles is less than 1." << endl;
		inputFile.close();
		return 0;
	}
	
	// allocation of arrays
	sconf.x = new float[sconf.amount];
	sconf.y = new float[sconf.amount];
	sconf.z = new float[sconf.amount];
	sconf.m = new float[sconf.amount];
	sconf.vx = new float[sconf.amount];
	sconf.vy = new float[sconf.amount];
	sconf.vz = new float[sconf.amount];
	
	for(int i = 0; i < sconf.amount; i++)
	{
		// read line
		inputLine = "";
		getline(inputFile,inputLine);
		
		#ifdef DEBUG_INPUT
		cout << "READ:" << inputLine << endl;
		#endif
		// create string stream from the line
		istringstream iss (inputLine);		
		
		#ifdef DEBUG_INPUT
		cout << "ISS CONTENT:" << iss.str() << endl;
		#endif
		// read values (floats)
		iss >> sconf.x[i] >> sconf.y[i] >> sconf.z[i];
		iss >> sconf.m[i];
		iss >> sconf.vx[i] >> sconf.vy[i] >> sconf.vz[i];
		
		// check eof - should be true with the last line read
		if(inputFile.eof())
		{
			#ifdef DEBUG_INPUT
			cout << "END OF FILE REACHED." << endl;
			#endif
			if((i+1) != sconf.amount) cout << "ERROR: Amount of particles in header is not equal to real amount of particles in the file." << endl;
			break;
		}
		
	}	
	
	// check output
	#ifdef DEBUG_INPUT
	for(int i = 0; i < amount; i++)
	{
		cout << "POS: [" << sconf.x[i] << ';' << sconf.y[i] << ';' << sconf.z[i] << "] | MASS: " << 
			sconf.m[i] << " | VELOCITIES: [" << sconf.vx[i] << ';' << sconf.vy[i] << ';' << sconf.vz[i] << ']' << endl;
	}
	#endif
	
	// close file
	inputFile.close();
	
	return 0;
}


/**
 * Create text file with particles.
 * Expects filled (already generated) arrays for x, y, ...
 */
int createOutputFile(const char * outputFilename, SimConfig & sconf)
{
	// open file
	ofstream outputFile(outputFilename, ios::out);
	
	// check if the file is opened
	if(!outputFile.is_open())
	{
		cout << "ERROR: Could not open the output file." << endl;
		return 0;
	}
	
	// expected format of file:
	/* (first line of file): amount (= integer, amount of particles stored in this file)
	 * (every other line): x y z m vx vy vz<lf> (x y (...) vz = floats)
	 */
	
	// first line = amount
	outputFile << sconf.amount << ' ' << sconf.simulation_steps << ' '
				<< sconf.width << ' '  << sconf.height << ' '  << sconf.depth << ' '
				 << sconf.max_speed << ' ' << sconf.max_weight << "\n";
	
	// one particle per one line
	for(int i = 0; i < sconf.amount; i++)
	{
		outputFile << sconf.x[i] << ' ' << sconf.y[i] << ' ' << sconf.z[i] << ' ';
		outputFile << sconf.m[i] << ' ';
		outputFile << sconf.vx[i] << ' ' << sconf.vy[i] << ' ' << sconf.vz[i] << '\n';
	}
	
	outputFile.close();	
	
	return 0;
}


#ifdef DEBUG_TEST_IO
int main(int argc, char **argv)
{
	SimConfig sconf;
	sconf.amount = 6;
	
	sconf.x = new float[6];
	sconf.y = new float[6];
	sconf.z = new float[6];
	
	sconf.m = new float[6];
	
	sconf.vx = new float[6];
	sconf.vy = new float[6];
	sconf.vz = new float[6];
	
	for(int i = 1; i < 7; i++)
	{
		sconf.x[i-1] = sconf.y[i-1] = sconf.z[i-1] = sconf.m[i-1] = i*10;
		sconf.vx[i-1] = sconf.vy[i-1] = sconf.vz[i-1] = i;
	}

	createOutputFile("testIO.txt", sconf);
	
	processInputFile("testIO.txt");
	
	return 0;
}
#endif

