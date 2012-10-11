#ifndef __FFTUTILS__
#define __FFTUTILS__

// Global Vars

//*******3dFFT configuration file setting STARTS************//
extern unsigned blockSize;
extern unsigned fftAlgo;
extern unsigned print;

extern unsigned xRange;
extern unsigned yRange;
extern unsigned zRange;

//*******3dFFT configuration file setting ENDS************//

// Global Vars end


// Functions Prototype 

void cleanup();
bool readConfig(const char* const fName);
unsigned initExecution(const unsigned size, const unsigned samplesize);
int allocateHostMemory(const unsigned size, const unsigned samplesize);
// Functions Prototype end

#endif
