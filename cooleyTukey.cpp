#include <iostream>
#include "cooleyTukey.h"
#include "fftutils.h"

using namespace std;


int initialxRange;
int initialyRange;
int initialzRange;

/* Functions Implementation*/

bool runCooleyTukey(const char* const argv[] ){
	
	const unsigned n = xRange;
	const unsigned size = xRange*yRange*zRange;
	int ASPAN = xRange*yRange*zRange;
	int zBar, yBar, xBar;
	int full3dfft=1;

	initialxRange = xRange;
	initialyRange = yRange;
	initialzRange = zRange;

	if (zRange > 1) {
		zBar = xRange; yBar = yRange; xBar = zRange;
	}
	else { 
		zBar = zRange; yBar = xRange; xBar = yRange;
	}
	
	int show_result = 0;
	int FFT_type = 0;
	
	printf("hello From runCooleyTukey\n");

	if (!initExecution(size, n)) {
		return false;
	}

	// matrix dimension
	unsigned sizeOnCPU = size;
	return true;

}

/* Functions Implementation end*/
