#include <stdio.h>
#include <fstream>
#include <iostream>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include "fftutils.h"

using namespace std;

// Global Vars

unsigned blockSize;
unsigned fftAlgo;
unsigned print;

//*******3dFFT default configuration file setting STARTS************//
unsigned xRange = 4;
unsigned yRange = 4;
unsigned zRange = 4;

unsigned xFRange = xRange/2;
unsigned yFRange = yRange/2;
unsigned zFRange = zRange/2;
//*******3dFFT default configuration file setting ENDS************//


// host memory
// h_Freal and h_Fimag represent the input signal to be transformed.
// h_Rreal and h_Rimag represent the transformed output.
double*  h_Freal = 0;
double*  h_Fimag = 0;
double*  h_Rreal = 0;
double*  h_Rimag = 0;

///////////////////// 3d fft host memory for N vector memeory allocation and initialization START /////////////////////////////////

double* hraVec = 0;
double* hiaVec = 0;

double* hrmVecI = 0;
double* hrmVecJ = 0;
double* hrmVecK = 0;
double* himVecI = 0;
double* himVecJ = 0;
double* himVecK = 0;

double* hrhVecI = 0;
double* hrhVecJ = 0;
double* hrhVecK = 0;
double* hihVecI = 0;
double* hihVecJ = 0;
double* hihVecK = 0;

//after bit reversal

double* hrRaVec = 0;
double* hiRaVec = 0;

double* hrRmVecI = 0;
double* hrRmVecJ = 0;
double* hrRmVecK = 0;
double* hiRmVecI = 0;
double* hiRmVecJ = 0;
double* hiRmVecK = 0;

double* hrRhVecI = 0;
double* hrRhVecJ = 0;
double* hrRhVecK = 0;
double* hiRhVecI = 0;
double* hiRhVecJ = 0;
double* hiRhVecK = 0;

///////////////////// 3d fft host memory for N vector memeory allocation and initialization START /////////////////////////////////

// Global Vars end

/* Functions Implementation*/

void cleanup(){
	cout<< "Cleanup"<<endl;

}


bool readConfig(const char* const fName){
	printf("Reading the file - %s - [start]\n",fName);
	ifstream file(fName);
    if (!file.is_open()) {
        cout << "Config file cannot be opened" << endl;
        return false;
    }
	cout << "Reading from config file....." << endl;
    while (!file.eof()) {
        string line;
        getline(file, line);
        if (line[0] == '\0') break;
        if (line[0] == '#') continue; 
        char config[13];
        unsigned val = 0;
        sscanf(line.c_str(), "%s %d", config, &val);
		if (!strcmp(config, "BLOCK_SIZE")) {
            if (val == 0) {
                cout << "Block size cannot be zero" << endl;
            }   
            blockSize = val;
        } else if (!strcmp(config, "FFT_ALGO")) {
            if (val > 5) {
                cout << "FFT_ALGO config should be from 1 to 5." << endl;
                return false;
            }
            fftAlgo = val;
		}else if (!strcmp(config, "xRange")) {
            if (val == 0) {
                cout << "xRange Should be greater than zero." << endl;
                return false;
            }
            xRange = val;
	 	}else if (!strcmp(config, "yRange")) {
            if (val == 0) {
                cout << "yRange Should be greater than zero." << endl;
                return false;
            }
            yRange = val;
	 	}else if (!strcmp(config, "zRange")) {
            if (val == 0) {
                cout << "zRange Should be greater than zero." << endl;
                return false;
            }
            zRange = val;	
	
        }  else if (!strcmp(config, "PRINT_RESULT")) {
            print = val;
        }   else {
            cout << "Invalid config " << config << endl;
            return false;
        } 
    }

    return true;
}

unsigned initExecution(const unsigned size, const unsigned samplesize){

    // Allocate host memory
    allocateHostMemory(size, samplesize);
    /*
    if (deviceCount) {
        cout << "Initializing device(s).." << endl;
        // create the OpenCL context on available GPU devices
        init_cl_context(CL_DEVICE_TYPE_GPU);

        const cl_uint ciDeviceCount = getDeviceCount();


        if (!ciDeviceCount) {
            printf("No opencl specific devices!\n");
            return 0;
        }

        printf("Creating Command Queue...\n");
        // create a command queue on device 1
        for (unsigned i = 0; i < deviceCount; ++i) {
            createCommandQueue(i);
        }
    }
    */
    return 1;
}

// HOST memory initialization

int allocateHostMemory(const unsigned size, const unsigned n){

    //memory alloc for h_Freal,h_Fimag,h_Rreal,h_Rimag,
    h_Freal = (double *) malloc(sizeof(double) * size);
    h_Fimag = (double *) malloc(sizeof(double) * size);
    h_Rreal = (double *) malloc(sizeof(double) * size);
    h_Rimag = (double *) malloc(sizeof(double) * size);
    
    //initialization for h_Freal,h_Fimag,h_Rreal,h_Rimag,
    for (unsigned i = 0 ; i < size; ++i) {

        h_Freal[i] = i+1.0;
        h_Fimag[i] = 0.0;

        h_Rreal[i] = 0.0;
        h_Rimag[i] = 0.0;
    } 

    ///////////////////// 3d fft host memory for N or aVec mVec and hVec vectors memeory allocation and initialization START /////////////////////////////////

    int xReal = 2;
    int yReal = 2;
    int zReal = 2;

    double MX = 26.726124191242441;
    double MY = 53.452248382484882; 
    double MZ = 80.178372573727316;

    printf("Real Dims: (%d,%d,%d)\n", xRange/2, yRange/2, zRange/2);
    //printf("FFTs Dims: (%d,%d,%d)\n", xFRange, yFRange, zFRange);
    printf("Finl Dims: (%d,%d,%d)\n", xRange, yRange, zRange);
    
    static const char *path[7];
    
    path[0] = "./datainput/1.log";
    path[1] = "./datainput/2.log"; 
    path[2] = "./datainput/3.log"; 
    path[3] = "./datainput/4.log"; 
    path[4] = "./datainput/5.log"; 
    path[5] = "./datainput/6.log";
    path[6] = "./datainput/224.omf"; 

    // [start] aVec matrix, one for each log file 
    
    //// before bit reversal
    hraVec = (double *) malloc(sizeof(double)*xRange*yRange*zRange*6);
    hiaVec = (double *) malloc(sizeof(double)*xRange*yRange*zRange*6);
    //// after bit reversal
    hrRaVec = (double *) malloc(sizeof(double)*xRange*yRange*zRange*6);
    hiRaVec = (double *) malloc(sizeof(double)*xRange*yRange*zRange*6);
    
    // [end] aVec matrix, one for each log file 

    // [start] mVec matrix of double

    //// before bit reversal
    hrmVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrmVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrmVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    himVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    himVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    himVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    //// after bit reversal
    hrRmVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrRmVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrRmVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hiRmVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hiRmVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hiRmVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);

    // [end] mVec matrix of double

    // [start] hVec matrix of double

    //// before bit reversal
    hrhVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrhVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrhVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hihVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hihVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hihVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    //// after bit reversal
    hrRhVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrRhVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hrRhVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hiRhVecI = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hiRhVecJ = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    hiRhVecK = (double *) malloc(sizeof(double)*xRange*yRange*zRange);
    // [end] hVec matrix of double

    // Read from log file 1.log 2.log 3.log 4.log 5.log 6.log
    char buffer[100];
    char *buffer_inter;
    char *ptr = buffer;
    char *remain;
    FILE *fd[7];
    int ASPAN = xRange*yRange*zRange;

    int x,y;
    for (x = 0; x < 6; x++) {
        fd[x] = fopen(path[x], "r");
        if (fd[x] != NULL) {
            y = 0;
            while( fgets(buffer, sizeof(buffer), fd[x] ) != NULL) 
            {
                hraVec[(x*ASPAN+y)] = atof(buffer);
                hiaVec[(x*ASPAN+y)] = 0.0;

                hrRaVec[(x*ASPAN+y)] = 0.0;
                hiRaVec[(x*ASPAN+y)] = 0.0;

                y++;
            }
        }
        else {
            printf("File Does not Exist\n");
            return 1;
        }
        //fclose(fd[0]); fiz
        fclose(fd[x]); //ki fix
    }

    
    // printf("aVEC STARTING VALUES \n");
    // printMe(0,hraVec, hiaVec, zRange, yRange, xRange, 0*ASPAN );

    //******************************** mVec Initialization  STARTS ***************************************

    int XSTART, ZSTART;
    int z;
    /*for (z = 0; z < zRange * xRange * yReal; z++) {
        mVecI[2*z] = cos((2*pi*z)/(2*xReal));
        mVecI[2*z+1] = -1*sin((2*pi*z)/(2*xReal));
    }*/

    fd[6] = fopen(path[6], "r");
    if (fd[6] != NULL ) {
        for (z = 0; zRange == 1 ? z < zRange : z < zRange/2; z++) {
            ZSTART = xRange * yRange * z;
            for (y = 0; y < yRange; y++) {
                XSTART = y * xRange;
                if (y > (yRange/2-1) && y > 1) {
                    for (x = 0; x < xRange; x++) {
                        hrmVecI[(ZSTART+XSTART+x)] = 0;
                        himVecI[(ZSTART+XSTART+x)] = 0;
                        hrmVecJ[(ZSTART+XSTART+x)] = 0;
                        himVecJ[(ZSTART+XSTART+x)] = 0;
                        hrmVecK[(ZSTART+XSTART+x)] = 0;
                        himVecK[(ZSTART+XSTART+x)] = 0;
                    }
                }
                else {
                    for (x = 0; x < xRange; x++) {
                        if (x > (xRange/2 - 1) && x > 1) {
                            hrmVecI[(ZSTART+XSTART+x)] = 0;
                            himVecI[(ZSTART+XSTART+x)] = 0;
                            hrmVecJ[(ZSTART+XSTART+x)] = 0;
                            himVecJ[(ZSTART+XSTART+x)] = 0;
                            hrmVecK[(ZSTART+XSTART+x)] = 0;
                            himVecK[(ZSTART+XSTART+x)] = 0;
                        }
                        else {
                            // From File
                            fgets(buffer, sizeof(buffer), fd[6] );
                            ptr = buffer;
                            buffer_inter = strtok_r(ptr, "\t", &remain);
                            hrmVecI[(ZSTART+XSTART+x)] = atof(ptr);//MX;
                            himVecI[(ZSTART+XSTART+x)] = 0;
                            ptr = remain;
                            
                            buffer_inter = strtok_r(ptr, "\t", &remain);
                            hrmVecJ[(ZSTART+XSTART+x)] = atof(ptr);//MY;
                            himVecJ[(ZSTART+XSTART+x)] = 0;
                            ptr = remain;

                            buffer_inter = strtok_r(ptr, "\n", &remain);
                            hrmVecK[(ZSTART+XSTART+x)] = atof(ptr);//MZ;
                            himVecK[(ZSTART+XSTART+x)] = 0;
                            ptr = remain;
                        }
                    }
                }
            }
        }
        if (z > 1) {
            for (z = zRange/2; z < zRange; z++) {
                ZSTART = xRange * yRange * z;
                for (y = 0; y < yRange; y++) {
                    XSTART = y * xRange;
                    if (y > (yRange/2-1) && y > 1) {
                        for (x = 0; x < xRange; x++) {
                            hrmVecI[(ZSTART+XSTART+x)] = 0;
                            himVecI[(ZSTART+XSTART+x)] = 0;
                            hrmVecJ[(ZSTART+XSTART+x)] = 0;
                            himVecJ[(ZSTART+XSTART+x)] = 0;
                            hrmVecK[(ZSTART+XSTART+x)] = 0;
                            himVecK[(ZSTART+XSTART+x)] = 0;
                        }
                    }
                    else {
                        for (x = 0; x < xRange; x++) {
                            if (x > (xRange/2 - 1) && x > 1) {
                                hrmVecI[(ZSTART+XSTART+x)] = 0;
                                himVecI[(ZSTART+XSTART+x)] = 0;
                                hrmVecJ[(ZSTART+XSTART+x)] = 0;
                                himVecJ[(ZSTART+XSTART+x)] = 0;
                                hrmVecK[(ZSTART+XSTART+x)] = 0;
                                himVecK[(ZSTART+XSTART+x)] = 0;
                            }
                            else {
                                hrmVecI[(ZSTART+XSTART+x)] = 0;
                                himVecI[(ZSTART+XSTART+x)] = 0;
                                hrmVecJ[(ZSTART+XSTART+x)] = 0;
                                himVecJ[(ZSTART+XSTART+x)] = 0;
                                hrmVecK[(ZSTART+XSTART+x)] = 0;
                                himVecK[(ZSTART+XSTART+x)] = 0;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(fd[6]);    
}
/* Functions Implementation end*/
