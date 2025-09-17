#include <cuda_runtime.h>
#include "../../lib/gpu_comm.cuh"

using std::cout ;


int main( int argc, char **argv ) {
  cudaDeviceProp GPU_PROPS ;

  printf( "=============================================================\n"
	  " Test script to get some features of the GPU.\n" ) ;

  GPU_query( &GPU_PROPS ) ;
  printf("=============================================================\n");
  exit(EXIT_SUCCESS);
}
