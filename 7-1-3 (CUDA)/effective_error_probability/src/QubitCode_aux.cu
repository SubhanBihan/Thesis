// Need to change encoded logical qubit (output of getLogicalKet) according to the QECC

#include <cmath>

#include "QubitCode.cuh"

std::vector<datatype> QubitCode::getLogicalKet(const double theta) {
    const datatype s = std::sin(0.5 * theta);
    const datatype c = std::cos(0.5 * theta);

    // 0*c + 1*s, where 0 and 1 are logical codewords for Stabilizer code!
    std::vector<datatype> output = {c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0)};
    
    return output;
}

std::vector<datatype> QubitCode::getLogicalDensityMatrix(const double theta) {
    std::vector<datatype> cArrCodeword(FastPauli::CArrLen<NUMQUBITS>);
    const std::vector<datatype> ketL = QubitCode::getLogicalKet(theta);
    FastPauli::ketToCArr(cArrCodeword, ketL, NUMQUBITS);
    return cArrCodeword;
}
