#include <stdio.h>
#include "fftutils.h"
#include <iostream>

using namespace std;

int main(int argc, char **argv){


	//read from configuration file
	const char* const fName = (argc == 2) ? argv[1] : 0;

	if (fName) {
         if (!readConfig(fName)) {
              cout << "Error in config file abort" << endl;  
              return 0;
         }
         
    } else {
         cout << "Config file not specified. Executing with default config" << endl;
    }

    //print configuration file
    cout << "BLOCK_SIZE = " << blockSize << endl;
    cout << "PRINT_RESULT " << print << endl;
    cout << "xRange = " << xRange << endl;
    cout << "yRange = " << yRange << endl;
    cout << "zRange = " << zRange << endl;

    //print selected algorithm
    if (fftAlgo==2){
    	cout << "FFT algorithm = COOLEY_TUKEY" << endl;
    }else{
    	cout << "Unknow FFT algorithm" << endl;
    }
}
