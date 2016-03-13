#include <iostream> // cout, etc.
#include <fstream> // file
#include <sstream> // stringstream

using namespace std;

#define DEBUG_INPUT // debug output (to cout) of input process function
#define DEBUG_TEST_IO // run test (compile with main)

/**
 * Process input file - read particles from file.
 * TODO: Add references to x, y, z, ... (?maybe a structure?) as parameters so the input can be actually used.
 */
int processInputFile(const char * inputFilename)
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
	/* (first line of file): amount (= integer, amount of particles stored in this file)
	 * (every other line): x y z m vx vy vz<lf> (x y (...) vz = doubles)
	 */
	
	// amount of particles stored in file
	int amount = 0;
	
	string inputLine = "";
	
	// read the header - amount of particles in file
	getline(inputFile,inputLine);
	
	// create stringstream and read the amount (?and other possible parameters?)
	istringstream issHeader(inputLine);
	issHeader >> amount;
	
	// check value of amount - should be more than 0
	if(amount < 1) 
	{
		cout << "ERROR: Amount of particles less than 1." << endl;
		inputFile.close();
		return 0;
	}
	
	// allocation of arrays
	double * x = new double[amount];
	double * y = new double[amount];
	double * z = new double[amount];
	double * m = new double[amount];
	double * vx = new double[amount];
	double * vy = new double[amount];
	double * vz = new double[amount];
	
	for(int i = 0; i < amount; i++)
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
		// read values (doubles)
		iss >> x[i] >> y[i] >> z[i];
		iss >> m[i];
		iss >> vx[i] >> vy[i] >> vz[i];
		
		// check eof - should be true with the last line read
		if(inputFile.eof())
		{
			#ifdef DEBUG_INPUT
			cout << "END OF FILE REACHED." << endl;
			#endif
			if((i+1) != amount) cout << "ERROR: Amount of particles in header is not equal to real amount of particles in the file." << endl;
			break;
		}
		
	}	
	
	// check output
	#ifdef DEBUG_INPUT
	for(int i = 0; i < amount; i++)
	{
		cout << "POS: [" << x[i] << ';' << y[i] << ';' << z[i] << "] | MASS: " << m[i] << " | VELOCITIES: [" << vx[i] << ';' << vy[i] << ';' << vz[i] << ']' << endl;
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
int createOutputFile(const char * outputFilename, double * x, double * y, double * z, double * m, double * vx, double * vy, double * vz, int amount)
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
	 * (every other line): x y z m vx vy vz<lf> (x y (...) vz = doubles)
	 */
	
	// first line = amount
	outputFile << amount << '\n';
	
	// one particle per one line
	for(int i = 0; i < amount; i++)
	{
		outputFile << x[i] << ' ' << y[i] << ' ' << z[i] << ' ';
		outputFile << m[i] << ' ';
		outputFile << vx[i] << ' ' << vy[i] << ' ' << vz[i] << '\n';
	}
	
	outputFile.close();	
	
	return 0;
}


#ifdef DEBUG_TEST_IO
int main(int argc, char **argv)
{
	int amount = 6;
	
	double x[] = {10,20,30,40,50,60};
	double y[] = {10,20,30,40,50,60};
	double z[] = {10,20,30,40,50,60};

	double m[] = {10,20,30,40,50,60};

	double vx[] = {10,10,10,10,10,10};
	double vy[] = {1,2,3,4,5,6};
	double vz[] = {1,2,3,4,5,6};
	createOutputFile("testIO.txt", x,y,z,m,vx,vy,vz,amount);
	
	processInputFile("testIO.txt");
	
	return 0;
}
#endif

