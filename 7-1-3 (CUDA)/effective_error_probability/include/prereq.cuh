// Simple header to store libraries

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <unordered_map>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>  // Fix the path
#include <qpp/qpp.hpp>
#include "params.cuh"

extern std::vector<std::string> corrections;
extern std::vector<std::string> corrections_erasure;
extern std::vector<std::string> corrections_flag;

extern __constant__ datatype I_d[], X_d[], Y_d[], Z_d[];

typedef Eigen::Matrix<datatype, DIM, DIM> Operator_t;

// INITIALIZATION must always be in row-major format
const Eigen::Matrix2f I_h = {1, 0, 0, 1};
const Eigen::Matrix2f X_h = {0, 1, 1, 0};
const Eigen::Matrix2f Z_h = {1, 0, 0, -1};
const Eigen::Matrix2f Y_h = {0, -1, 1, 0};

const std::unordered_map<char, Eigen::Matrix2f> pauli_mat = {
    {'I', I_h},
    {'X', X_h},
    {'Y', Y_h},
    {'Z', Z_h}
};

extern std::vector<Operator_t> Stabilizers;

#define CHECK_CUDA(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error in %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

#define CHECK_CUBLAS(call) { \
    cublasStatus_t status = call; \
    if (status != CUBLAS_STATUS_SUCCESS) { \
        fprintf(stderr, "cuBLAS error in %s:%d: %d\n", __FILE__, __LINE__, status); \
        exit(EXIT_FAILURE); \
    } \
}

class cuda_global {
    public:
        std::shared_ptr<datatype> arr;  // Pointer to the device array. Switch to manual memory management later if needed
        size_t size;

        explicit cuda_global(size_t size);

        cuda_global(size_t size, datatype *host_arr);
        cuda_global(Eigen::Matrix<datatype, Eigen::Dynamic, Eigen::Dynamic> &host_arr);

        // These are only allowed because of shared_ptr
        // Must define ourselves if using manual memory management
        cuda_global(const cuda_global &other) = default;
        cuda_global(cuda_global &&other) = default;
        cuda_global &operator=(const cuda_global &other) = default;
        cuda_global &operator=(cuda_global &&other) = default;      
        ~cuda_global() = default;

        // Can be expanded further (operations etc.)
        
};

std::ostream &operator<<(std::ostream &os, const cuda_global &d_arr);

void copyUpperTriangularColumnMajor(const Eigen::Matrix<datatype, dynamic, dynamic>& mat, datatype* dest, unsigned reduced_size);
