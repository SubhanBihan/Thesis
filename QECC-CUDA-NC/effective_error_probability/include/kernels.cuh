#pragma once
#include "prereq.cuh"

__global__ void kron_product(datatype* A, int m, int n, datatype* B, int p, int q, datatype* C);

// __global__ void combine_and_scale_2d(datatype *C, datatype *X, datatype *Y, datatype *Z,
//                                      datatype alpha, datatype beta, unsigned width);

__global__ void combine_2d(datatype *C, datatype *X, datatype *Y, unsigned width);

__global__ void combine_and_scale(datatype *C, datatype **inputs, unsigned num_inputs, 
                                        datatype alpha, datatype beta, unsigned width);

__global__ void X_dm(datatype* output, datatype* input, unsigned k, unsigned N);

__global__ void Y_dm(datatype* output, datatype* input, unsigned k, unsigned N);

__global__ void Z_dm(datatype* output, datatype* input, unsigned k, unsigned N);

__global__ void CX_dm(datatype* output, datatype* input, 
                      unsigned control_qubit, unsigned target_qubit, unsigned N);

__global__ void CY_dm(datatype* output, datatype* input, 
                      unsigned control_qubit, unsigned target_qubit, unsigned N);
                      
__global__ void CZ_dm(datatype* output, datatype* input, 
                      unsigned control_qubit, unsigned target_qubit, unsigned N);

__global__ void trace_dm(datatype* output, datatype* input, unsigned N);

__global__ void fidel_dm(datatype* output, datatype* input, datatype* logDM, unsigned N);

namespace NoisyPauli {
    void NCG_dm(char G, datatype* output, datatype* input, const unsigned N, const unsigned controlQubit, const unsigned targetQubit, const datatype err,
                        std::vector<cudaStream_t>& streams, datatype** d_temps);
    // void NCX_dm(cuda_global& output, const cuda_global& input, const unsigned N, const unsigned controlQubit, const unsigned targetQubit, const datatype err,
    //             std::vector<cudaStream_t>& streams, datatype** d_temps);
    // void NCZ_dm(cuda_global& output, const cuda_global& input, const unsigned N, const unsigned controlQubit, const unsigned targetQubit, const datatype err,
    //             std::vector<cudaStream_t>& streams, datatype** d_temps);
};