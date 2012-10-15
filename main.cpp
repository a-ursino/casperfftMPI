#include <stdio.h>
#include <iostream>
#include "fftutils.h"
#include "cooleyTukey.h"
#include "mpi.h"

using namespace std;

int root=0;
int initialxRange;
int initialyRange;
int initialzRange;


int main(int argc, char **argv){

  bool result = true;

  MPI_Init(NULL,NULL);

  //get the number of processes
  int world_size=0;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);

  //get the rank
  int rankid=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rankid);
	// Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  if(rankid==root){
    printf("Hello world from processor %s, rank %d out of %d processors. THE MASTER\n", processor_name, rankid, world_size);

    //read from configuration file
  	const char* const fName = (argc == 2) ? argv[1] : 0;

  	if (fName) {
      if (!readConfig(fName)) {
        printf("Error in config file abort");
        return 0;
      }else{
        printf("Reading the file - %s - [end]\n", fName);        
      }
           
    } else {
      printf("Config file not specified. Executing with default config");
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
    }else{
      printf("Hello world from processor %s, rank %d out of %d processors. SLAVE\n", processor_name, rankid, world_size);

    }

    //set up enviroment....
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
    unsigned sizeOnCPU=size; // matrix dimension

    //...and allocate memory on root node
    if (rankid==root){
      if (!initExecution(size, n)) {
        return false;
      }      
    }

    if(full3dfft==1){




    }

    //All

    MPI_Finalize();
}
