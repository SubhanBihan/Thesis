// Simple header to store libraries

#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <iterator>
#include <unordered_map>
#include <cuda_runtime.h>
#include <cublas_v2.h>
// #include <Eigen/Dense>
// #include <unsupported/Eigen/KroneckerProduct>
// #include <qpp/qpp.hpp>
#include "params.cuh"

// extern __constant__ datatype DM_H1P[4], DM_H2P0[16];
// extern std::vector<datatype*> DM_H1;
// extern datatype* DM_H2P0_ptr;

void readLinesFromFile(std::vector<std::string>& lines, const std::string &filename);

// typedef Eigen::Matrix<datatype, DIM, DIM> Operator_t;
// typedef Eigen::Matrix<datatype, Eigen::Dynamic, Eigen::Dynamic> Dynamic_t;

extern std::vector<std::string> corrections;
extern std::vector<std::string> corrections_erasure;
extern std::vector<std::string> corrections_flag;
extern std::vector<std::string> stabilizers;

// extern std::vector<Operator_t> Stabilizers;

#define CHECK_CUDA(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error in %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

#define CUDA_KERNEL_CHECK() \
    do { \
        cudaError_t err = cudaGetLastError(); \
        if (err != cudaSuccess) { \
            fprintf(stderr, "Kernel launch error %s at %s:%d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

#define CHECK_CUBLAS(call) \
    do { \
        cublasStatus_t status = call; \
        if (status != CUBLAS_STATUS_SUCCESS) { \
            fprintf(stderr, "cuBLAS error in %s:%d: %d\n", __FILE__, __LINE__, status); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

void SyncStreams(std::vector<cudaStream_t>& streams);

class cuda_global {
    private:
        datatype *arr;  // Pointer to the device array.

    public:
        size_t size;

        cuda_global(){};
        explicit cuda_global(size_t size);

        // cuda_global(size_t size, datatype *host_arr);
        cuda_global(const std::vector<datatype>& host_vec);

        cuda_global(const cuda_global &other) = delete;
        cuda_global(cuda_global &&other);
        cuda_global &operator=(const cuda_global &other);
        cuda_global &operator=(cuda_global &&other);      
        ~cuda_global();

        datatype *get(void) const;
        void SetZero(void);
        friend std::ostream &operator<<(std::ostream &os, const cuda_global &d_arr);
        friend void swap(cuda_global &first, cuda_global &second);
        void copy_to_host(std::vector<datatype>& host_vec) const;

};

extern dim3 blockSize;
dim3 getGridSize(const unsigned mat_rows, const unsigned mat_cols);

extern std::vector<cuda_global> ancilla, ancilla_flag;     // Too big for constant memory

extern cuda_global DM_H1P, DM_H2P0;