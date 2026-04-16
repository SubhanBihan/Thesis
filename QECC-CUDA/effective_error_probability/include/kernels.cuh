#pragma once
#include "prereq.cuh"

__global__ void kron_product(datatype* A, int m, int n, datatype* B, int p, int q, datatype* C);

// __global__ void combine_and_scale_2d(datatype *C, datatype *X, datatype *Y, datatype *Z,
//                                      datatype alpha, datatype beta, unsigned width);

__global__ void combine_2d(datatype *C, datatype *X, datatype *Y, unsigned width);

__global__ void combine_and_scale(datatype *C, datatype **inputs, unsigned num_inputs, 
                                        datatype alpha, datatype beta, unsigned width);

__global__ void combine_and_scale_3(datatype *C, datatype *X, datatype *Y, datatype *Z, 
                                    datatype alpha, datatype beta, unsigned width);

__global__ void combine_and_scale_12(datatype *C, datatype *d1, datatype *d2, datatype *d3, datatype *d4, datatype *d5, datatype *d6, datatype *d7, datatype *d8, datatype *d9, datatype *d10, datatype *d11, datatype *d12, 
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
    void NCG_dm(char G, cuda_global& output, cuda_global& input, const unsigned N, const unsigned controlQubit, const unsigned targetQubit, const datatype err,
                std::vector<cudaStream_t>& streams, std::vector<cuda_global>& d_temps);
};