#include <stdio.h>
#include <iostream>
#include <sstream>
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

    //set up enviroment and initialize shared variable

    const unsigned n = xRange;
    const unsigned size = xRange*yRange*zRange;
    int ASPAN = xRange*yRange*zRange;
    int zBar, yBar, xBar;
    int full3dfft=1;
    int show_result = 0;
    int FFT_type = 0;
    unsigned matrix_size=size; // matrix dimension sizeOnCPU
    int num_elements_per_proc=size; //matrix dimension

    initialxRange = xRange;
    initialyRange = yRange;
    initialzRange = zRange;

    if (zRange > 1) {
      zBar = xRange; yBar = yRange; xBar = zRange;
    }
    else { 
      zBar = zRange; yBar = xRange; xBar = yRange;
    }



    //...and allocate memory on root node
    if (rankid==root){
      if (!initExecution(size, n)) {
        return false;
      }      
    }
    //buffer to receive data from scatter
    double *recv_hraVec_buffer=(double*)malloc(sizeof(double)*num_elements_per_proc);
    double *recv_hiaVec_buffer=(double*)malloc(sizeof(double)*num_elements_per_proc);

    if(full3dfft==1){
      show_result = 0;
      FFT_type = 0; //Forward FFT
      //----------------------------(2. aVec XX,XY,XZ,YY,YZ,ZZ [start])----------------------------------

      MPI_Scatter(hraVec, num_elements_per_proc, MPI_DOUBLE, recv_hraVec_buffer, num_elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Scatter(hiaVec, num_elements_per_proc, MPI_DOUBLE, recv_hiaVec_buffer, num_elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      //stringstream msg;
      //msg << "recv_hraVec_buffer recv_hiaVec_buffer from rankid " << rankid;
      //printMeInfo(msg.str(),0,recv_hraVec_buffer, recv_hiaVec_buffer, zRange, yRange, xRange, 0*ASPAN );

      //local buffer
      double *local_hrRaVec_buffer=(double*)malloc(sizeof(double)*num_elements_per_proc);
      double *local_hiRaVec_buffer=(double*)malloc(sizeof(double)*num_elements_per_proc);

      cooleyTukeyCpu3DFFT(0, n, matrix_size,recv_hraVec_buffer,recv_hiaVec_buffer,local_hrRaVec_buffer,local_hiRaVec_buffer,0,show_result,FFT_type,xRange,yRange,zRange);

      MPI_Gather(local_hrRaVec_buffer, num_elements_per_proc, MPI_DOUBLE, hrRaVec, num_elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gather(local_hiRaVec_buffer, num_elements_per_proc, MPI_DOUBLE, hiRaVec, num_elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);


      //----------------------------(2. aVec XX,XY,XZ,YY,YZ,ZZ [Ends])----------------------------------
      //all the slave must have complete

      //----------------------------(3 mVecI [Starts])----------------------------------
      if (rankid==root){
        int dest=1;
        printf("Sending the mVecI data from root to dest: %d ...\n", dest);
        MPI_Send(hrmVecI, matrix_size, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
        MPI_Send(himVecI, matrix_size, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
        printf("Sent the mVecI data from root to dest: %d end\n", dest);
      }
      if (rankid==1){
        int src=0;
        double *recv_hrmVecI_buffer=(double*)malloc(sizeof(double)*matrix_size);
        double *recv_himVecI_buffer=(double*)malloc(sizeof(double)*matrix_size);
        
        double *local_hrRmVecI_buffer=(double*)malloc(sizeof(double)*num_elements_per_proc);
        double *local_hiRmVecI_buffer=(double*)malloc(sizeof(double)*num_elements_per_proc);
        
        printf("Process 1: Receiving the mVecI data from root...\n");
        MPI_Recv(recv_hrmVecI_buffer, matrix_size, MPI_DOUBLE, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(recv_himVecI_buffer, matrix_size, MPI_DOUBLE, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 1: Received the mVecI data from root [end]\n");

        printf("Process 1: Processing with ct3dfft mVecI data ...\n");
        cooleyTukeyCpu3DFFT(0, n, matrix_size,recv_hrmVecI_buffer,recv_himVecI_buffer,local_hrRmVecI_buffer,local_hiRmVecI_buffer,0,show_result,FFT_type,xRange,yRange,zRange);
        printf("Process 1: Processed with ct3dfft mVecI data \n");
        //send the data processed back to root   
        printf("Sending the mVecI data processed from process 1 to root...\n");     
        MPI_Send(local_hrRmVecI_buffer, matrix_size, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
        MPI_Send(local_hiRmVecI_buffer, matrix_size, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
        printf("Sent the mVecI data processed from process 1 to root end\n");
      }
      //----------------------------(3 mVecI [Ends])----------------------------------

      //----------------------------(4 mVecJ [Starts])----------------------------------


      //----------------------------(4 mVecJ [Ends])----------------------------------

      //----------------------------(5 mVecK [Starts])----------------------------------


      //----------------------------(5 mVecK [Ends])----------------------------------

      if (rankid==root){
        //mVecI
        MPI_Recv(hrRmVecI, matrix_size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(hiRmVecI, matrix_size, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //mVecJ

        //mVecK
      }


      MPI_Barrier(MPI_COMM_WORLD);


      //convolveCPU(offset1,ASPAN);

    }

    //All

    MPI_Finalize();
}
