// Simple header to store number of qubits (n), 2^n, and other values
// Will be used in all other files that need these values

#ifndef NUMQUBITS

#include <vector>
#include <string>
#include <cmath>

#define NUMQUBITS 7
#define DIM 128
#define NUM_ANCILLA 6
#define DIM_ANCILLA 64

#define NUMPOINTS 10000
#define MAX_DIST 10000      // In km
// #define NUMTHREADS 16

using dtype = double;

extern std::vector<std::string> corrections;
extern std::vector<std::string> corrections_erasure;
extern std::vector<std::string> corrections_flag;

extern const dtype num;

#endif