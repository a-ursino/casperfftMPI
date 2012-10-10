#include <stdio.h>
#include "fftutils.h"
#include <iostream>

using namespace std;


unsigned blockSize;
unsigned fftAlgo;
unsigned print;

unsigned xRange;
unsigned yRange;
unsigned zRange;


/* Functions Implementation*/

void cleanup(){
	cout<< "Cleanup"<<endl;

}


bool readConfig(const char* const fName){
	printf("Reading the file %s\n",fName);
	return true;
}

/* Functions Implementation*/
