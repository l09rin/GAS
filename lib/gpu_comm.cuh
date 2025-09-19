/**
 * @file gpu_comm.cuh
 * @brief Utility functions and macros for CUDA error checking and GPU queries.
 */

#pragma once
#include <cuda_runtime.h>
#define CUDA_RUNTIME_H

#include <iostream>
using std::cout ;
using std::endl ;

// Error-checking macro
#define CUDA_CHECK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess) {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        exit(code);
    }
}

// Declarations
int getCudaCoresPerSM(int major, int minor) ;
void GPU_query( cudaDeviceProp *props ) ;
template<typename Func> void findOptimalGrid( dim3 &grid , dim3 &block , Func kernel , int dataSize , int SmemSize ) ;

/****************************************************************************************/
/******************************        Definitions        *******************************/

/**
 * @brief Get number of CUDA cores per streaming multiprocessor (SM).
 *
 * @param major Major compute capability.
 * @param minor Minor compute capability.
 * @return Number of cores per SM.
 */
int getCudaCoresPerSM(int major, int minor) {
    // Defines cores per SM (Streaming Multiprocessor) for various compute capabilities
    // Note: This table needs to be updated as new architectures are released
    // Source: NVIDIA architecture whitepapers / CUDA programming guides

    if (major == 9 && minor == 0) {  // Hopper (H100, GH100)
      return 128;
    } else if (major == 8 && minor == 9) {         // Ada Lovelace (SM 8.9, RTX 40-series)
        return 128;
    } else if (major == 8 && minor == 6) {  // Ampere (GA10x, RTX 30-series consumer)
        return 128;
    } else if (major == 8 && minor == 0) {  // Ampere (GA100, A100)
        return 64;
    } else if (major == 7 && minor == 5) {  // Turing (RTX 20-series, GTX 16-series)
        return 64;
    } else if (major == 7 && (minor == 0 || minor == 2)) {  // Volta GV100 (7.0) / Turing TU102 datacenter (7.2)
        return 64;
    } else if (major == 6 && (minor == 0 || minor == 1)) {  // Pascal GP100 (6.0), GP10x (6.1)
        return 64;
    } else if (major == 6 && minor == 2) {  // Pascal GP10B (Tegra X2)
        return 128;
    } else if (major == 5) { // Maxwell
        return 128;
    } else if (major == 3) { // Kepler
        return 192;
    } else if (major == 2 && (minor == 1 || minor == 0)) { // Fermi
        return (minor == 1) ? 48 : 32;
    } else {
        std::cerr << "Unknown device type" << std::endl;
        return 0;
    }
}
/****************************************************************************************/

/**
 * @brief Query GPU properties.
 *
 * Prints device information to stdout.
 *
 * @param props Pointer to cudaDeviceProp structure to fill.
 */
void GPU_query( cudaDeviceProp *props ) {
  int deviceId ;

  cudaGetDevice( &deviceId ) ;
  cudaGetDeviceProperties( props , deviceId ) ;

  double totGmem = (double)(props->totalGlobalMem) ;
  double totL2cache = (double)(props->l2CacheSize) ;
  double totCmem = (double)(props->totalConstMem) ;
  double totSmem = (double)(props->sharedMemPerBlock) ;
  char u_Gmem[] = " b" ;
  char u_Cmem[] = " b" ;
  char u_Smem[] = " b" ;
  char u_L2c[] = " b" ;
  if( totGmem > 1024 ) totGmem /= 1024 , strcpy( u_Gmem , "Kb" ) ;
  if( totGmem > 1024 ) totGmem /= 1024 , strcpy( u_Gmem , "Mb" ) ;
  if( totCmem > 1024 ) totCmem /= 1024 , strcpy( u_Cmem , "Kb" ) ;
  if( totCmem > 1024 ) totCmem /= 1024 , strcpy( u_Cmem , "Mb" ) ;
  if( totSmem > 1024 ) totSmem /= 1024 , strcpy( u_Smem , "Kb" ) ;
  if( totSmem > 1024 ) totSmem /= 1024 , strcpy( u_Smem , "Mb" ) ;
  if( totL2cache > 1024 ) totL2cache /= 1024 , strcpy( u_L2c , "Kb" ) ;
  if( totL2cache > 1024 ) totL2cache /= 1024 , strcpy( u_L2c , "Mb" ) ;

  cout << "============= GPU resources =============" << endl ;
  cout << "Device Name : " << props->name << endl ;
  if( props->integrated ) cout << "Integrated." << endl ;
  else cout << "Dedicated." << endl ;
  cout << "Total Global Memory ("<<u_Gmem<<") / Clock Rate (MHz) : " << totGmem << " / " << (double)(props->clockRate)/1000 << endl ;
  cout << "Compute Capability : " << props->major << "." << props->minor << endl ;
  cout << "Compute Mode : " << props->computeMode << endl ;
  cout << "            -----------------            " << endl ;
  cout << "               Processors" << endl ;
  cout << "MultiProcessors : " << props->multiProcessorCount << endl ;
  cout << "Cuda cores : " << getCudaCoresPerSM(props->major, props->minor) * props->multiProcessorCount << endl ;
  cout << "Max Threads per Block : " << props->maxThreadsPerBlock << endl ;
  cout << "Max Block : " << props->maxThreadsDim[0] << " x " << props->maxThreadsDim[1] << " x " << props->maxThreadsDim[2] << endl ;
  cout << "Max Grid : " << props->maxGridSize[0] << " x " << props->maxGridSize[1] << " x " << props->maxGridSize[2] << endl ;
  cout << "Warp Size : " << props->warpSize << endl ;
  cout << "            -----------------            " << endl ;
  cout << "                  Memory" << endl ;
  cout << "Clock (MHz) : " << props->memoryClockRate/1000 << endl ;
  cout << "Bus width (bit) : " << props->memoryBusWidth << endl ;
  cout << "L2 ("<<u_L2c<<") : " << totL2cache << endl ;
  cout << "ConstMem ("<<u_Cmem<<") : " << totCmem << endl ;
  cout << "SharedMem ("<<u_Smem<<") per Block : " << totSmem << endl ;
  cout << "            -----------------            " << endl ;
  cout << "Overlap : " << props->deviceOverlap << " " << props->asyncEngineCount << " " << props->concurrentKernels << endl ;
  cout << "=========================================" << endl ;
}
/****************************************************************************************/

/**
 * @brief Find optimal CUDA kernel launch configuration.
 *
 * @tparam Func Callable kernel function.
 * @param func Kernel function to evaluate.
 * @param blockSize Optimal block size (output).
 * @param minGridSize Minimum grid size (output).
 * @param dynamicSMemSize Dynamic shared memory size in bytes.
 */
template<typename Func>
void findOptimalGrid( dim3 &grid , dim3 &block , Func kernel , int dataSize , int SmemSize ) {
  int blockSize = 0;   // The optimal block size
  int minGridSize = 0; // Minimum grid size needed to achieve the maximum occupancy
  // Get the optimal block size for VV Kernel (before last argument is for dynamic shared memory...)
  cudaOccupancyMaxPotentialBlockSize( &minGridSize , &blockSize ,
				      kernel , SmemSize , dataSize ) ;
  // Calculate the grid size needed to cover all elements
  int gridSize = ( dataSize + blockSize - 1 ) / blockSize ;
  block = dim3(blockSize) ;
  grid = dim3(gridSize) ;
  // grid = dim3(max(minGridSize, gridSize)) ;
  cout << "Recommended min Grid size: " << minGridSize << endl ;
  cout << "Actual Grid size used: " << grid.x << endl ;
  cout << "Block size used: " << block.x << endl ;
}
/****************************************************************************************/
