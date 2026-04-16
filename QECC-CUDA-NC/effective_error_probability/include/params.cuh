// Simple header to store parameters

#pragma once

#include <cmath>

#define NUMQUBITS 7
#define DIM 128
#define NUM_ANCILLA 6
#define DIM_ANCILLA 64
#define STABILIZER_GATES (NUM_ANCILLA - 2)

#define NUMPOINTS 10000
#define MAX_DIST 10000      // In km
// #define NUMTHREADS 16
#define BLOCK_DIM 32
#define DM_SIZE (DIM * DIM)

using datatype = float;

extern const datatype num;

extern datatype alpha, beta;

#define traceDM(handle, N, A, I, results_d, i) CHECK_CUBLAS(cublasSdot(handle, N, A, N + 1, I, N + 1, results_d + i))
#define fidelDM(handle, N, A, logDM, results_d, i) CHECK_CUBLAS(cublasSdot(handle, N * N, A, 1, logDM, 1, results_d + i))