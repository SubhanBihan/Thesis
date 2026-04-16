#include "cuda_runtime.h"
#include <cstdint>
#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

std::vector<int> AP_host(const int start, const int end, const int d)
{
    assert(d > 0);
    const int n = (end - start) / d + 1; // equivalent to std::ceil((end-start)/d).
    std::vector<int> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = start + i * d;
    }
    return out;
}

__global__ void AP_dev(const int start, const int n, const int d, int* out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < n) out[i] = start + i * d; 
}

double get_time() {
    #ifdef _WIN32
        // Windows implementation
        LARGE_INTEGER frequency, counter;
        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&counter);
        return (double)counter.QuadPart / (double)frequency.QuadPart;
    #else
        // Unix/Linux/Mac implementation
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    #endif
}

int main()
{
    const int start = 0, end = 10000, d = 1;

    const int n = (end - start) / d + 1;

    int* AP_dev_out;
    std::vector<int> AP_dev_out_host(n);
    std::vector<int> AP_host_out(n);

    const int blockSize = 128;
    const int gridSize = ceil(n / blockSize);

    
    cudaMalloc(&AP_dev_out, n * sizeof(int));

    // Benchmark CPU implementation
    printf("Benchmarking CPU implementation...\n");
    double cpu_total_time = 0.0;
    for (int i = 0; i < 20; i++) {
        double start_time = get_time();

        AP_host_out = AP_host(start, end, d);

        double end_time = get_time();
        cpu_total_time += end_time - start_time;
    }
    double cpu_avg_time = cpu_total_time / 20.0;

    printf("Benchmarking GPU implementation...\n");
    double gpu_total_time = 0.0;
    double gpu_comp_time = 0.0;
    for (int i = 0; i < 20; i++) {
        double start_time = get_time();

        double gpu_comp_start = get_time();
        AP_dev<<<gridSize, blockSize>>>(start, n, d, AP_dev_out);
        double gpu_comp_end = get_time();
        cudaMemcpy(AP_dev_out_host.data(), AP_dev_out, n * sizeof(int), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        double end_time = get_time();
        gpu_total_time += end_time - start_time;
        gpu_comp_time += gpu_comp_end - gpu_comp_start;
    }
    double gpu_avg_time = gpu_total_time / 20.0;
    double gpu_avg_comp_time = gpu_comp_time / 20.0;
    
    // Print results
    printf("CPU average time: %f milliseconds\n", cpu_avg_time*1000);
    printf("GPU average time: %f milliseconds\n", gpu_avg_time*1000);
    printf("GPU computation average time: %f milliseconds\n", gpu_avg_comp_time*1000);
    printf("Speedup: %fx\n", cpu_avg_time / gpu_avg_time);

    cudaFree(AP_dev_out);
    return 0;
}

