#ifndef __FFTUTILS__
#define __FFTUTILS__

// Global Vars
extern unsigned blockSize;
extern unsigned fftAlgo;
extern unsigned print;
extern unsigned _show_result;

//*******3dFFT configuration file setting STARTS************//
extern unsigned xRange;
extern unsigned yRange;
extern unsigned zRange;

extern unsigned xFRange;
extern unsigned yFRange;
extern unsigned zFRange;
//*******3dFFT configuration file setting ENDS************//

extern double*  h_Freal;
extern double*  h_Fimag;
extern double*  h_Rreal;
extern double*  h_Rimag;

///////////////////// 3d fft host memory for N vector memeory allocation and initialization START /////////////////////////////////

extern double *hraVec;
extern double *hiaVec;

extern double *hrmVecI;
extern double *hrmVecJ;
extern double *hrmVecK;
extern double *himVecI;
extern double *himVecJ;
extern double *himVecK;

extern double *hrhVecI;
extern double *hrhVecJ;
extern double *hrhVecK;
extern double *hihVecI;
extern double *hihVecJ;
extern double *hihVecK;

/////after bit reversal

extern double *hrRaVec;
extern double *hiRaVec;

extern double *hrRmVecI;
extern double *hrRmVecJ;
extern double *hrRmVecK;
extern double *hiRmVecI;
extern double *hiRmVecJ;
extern double *hiRmVecK;

extern double *hrRhVecI;
extern double *hrRhVecJ;
extern double *hrRhVecK;
extern double *hiRhVecI;
extern double *hiRhVecJ;
extern double *hiRhVecK;

///////////////////// 3d fft host memory for N vector memeory allocation and initialization START /////////////////////////////////

// Global Vars end


// Functions Prototype 

void cleanup();
bool readConfig(const char* const fName);
unsigned initExecution(const unsigned size, const unsigned samplesize);
int allocateHostMemory(const unsigned size, const unsigned samplesize);

void printMe(const unsigned offset, double *datar,double *datai, int zR, int yR, int xR, int off);
void printMeInfo(const std::string& msg,const unsigned offset, double *datar,double *datai, int zR, int yR, int xR, int off);


// Functions Prototype end

#endif
