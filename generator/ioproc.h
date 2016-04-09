#include "SimConfig.h"

int processInputFile(const char * inputFilename, SimConfig & sconf);

int createOutputFile(const char * outputFilename, SimConfig & sconf);

int createGnuplotFile(const char * outputFilename, const char * outputDataFilename, SimConfig & sconf);
