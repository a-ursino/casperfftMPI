#include <stdio.h>
#include <iostream>
#include "fftutils.h"
#include "cooleyTukey.h"
#include "mpi.h"

using namespace std;

int root=0;


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
        result = runCooleyTukey(argv);
      }else{
      	cout << "Unknow FFT algorithm" << endl;
      }
    }else{
      printf("Hello world from processor %s, rank %d out of %d processors. SLAVE\n", processor_name, rankid, world_size);

    }

    //All

    float *rand_nums, *sub_rand_nums,*recv_buffer;
    int num_proc=6;
    int num_elements_per_proc=2;
    if (rankid==root){
      rand_nums= (float*)malloc(sizeof(float)*num_proc*num_elements_per_proc);
      for (int i = 0; i < num_elements_per_proc*num_proc; ++i){
        rand_nums[i]=2.5;
      }
    }
    
    sub_rand_nums=(float*)malloc(sizeof(float)*num_elements_per_proc);
    //MPI_Scatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator)
    MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums, num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    //cooleyTukeyCpu3DFFT(start, n, sizeOnCPU,hraVec,hiaVec,hrRaVec,hiRaVec,ASPAN,show_result,FFT_type);

    for (int i = 0; i < num_elements_per_proc; ++i){
      sub_rand_nums[i]=rankid;
      //printf("Processor %s, rank %d out of %d processors. [%f ]\n",processor_name, rankid, world_size,sub_rand_nums[i] );
    }

    if (rankid==root)
    {
      recv_buffer=(float*)malloc(sizeof(float)*num_proc*num_elements_per_proc);
    }
    //MPI_Gather(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator)
    MPI_Gather(sub_rand_nums, num_elements_per_proc, MPI_FLOAT, recv_buffer, num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rankid==root){   
      for (int i = 0; i < num_elements_per_proc*num_proc; ++i){      
        printf("Processor %s, rank %d out of %d processors. [%f ]\n",processor_name, rankid, world_size,sub_rand_nums[i] );
      }
    }





    MPI_Finalize();
}
