#include <cuda_runtime.h>
#include <iostream>


int main() {
    cudaDeviceProp prop;
    int device;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);

    std::cout << "Compute Capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "Number of Streaming Multiprocessors: " << prop.multiProcessorCount << std::endl;
    std::cout << "Max Threads per SM: " << prop.maxThreadsPerMultiProcessor << std::endl;
    std::cout << "Max Blocks per SM: " << prop.maxBlocksPerMultiProcessor << std::endl;
    std::cout << "Shared Memory per SM: " << prop.sharedMemPerMultiprocessor << " bytes" << std::endl;
    std::cout << "Registers per SM: " << prop.regsPerMultiprocessor << std::endl;
    std::cout << "Registers per Block: " << prop.regsPerBlock << std::endl;
    std::cout << "Total constant memory: " << prop.totalConstMem / 1024 << " KB" << std::endl;

    // L2 Cache size
    std::cout << "L2 Cache Size: " << prop.l2CacheSize << " bytes" << std::endl;


    int dev_count;  // Just use int directly since cudaGetDeviceCount expects int*
    cudaGetDeviceCount(&dev_count);

    std::cout << "Number of CUDA devices: " << dev_count << std::endl;

    std::cout << "Max threads per block along x: " << prop.maxThreadsDim[0] << std::endl;
    std::cout << "Warp size: " << prop.warpSize << std::endl;
    return 0;
}
