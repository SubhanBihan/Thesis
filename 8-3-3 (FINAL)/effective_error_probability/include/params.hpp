// Simple header to store number of qubits (n), 2^n, and other values
// Will be used in all other files that need these values

#ifndef NUMQUBITS

#include <vector>
#include <string>
#include <cmath>

#define NUMQUBITS 8
#define K 3
#define DIM 256
#define NUM_ANCILLA 5
#define DIM_ANCILLA 32

#define NUMPOINTS 100
#define MAX_DIST 10000      // In km
// #define NUMTHREADS 16

using dtype = double;

extern std::vector<std::string> corrections;
extern std::vector<std::string> corrections_erasure;
extern std::vector<std::string> corrections_flag;

extern const dtype num;

#endif