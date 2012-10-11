#include <stdio.h>
#include <fstream>
#include <iostream>
#include "fftutils.h"

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

/* Functions Implementation end*/
