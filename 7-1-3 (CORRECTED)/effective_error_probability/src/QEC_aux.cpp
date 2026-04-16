// Need to change encoded logical qubit (output of getLogicalKet)
// Changes have been made

#include "QEC.hpp"

std::vector<dtype> QEC::getLogicalKet(const dtype theta) {
    const dtype s = std::sin(0.5 * theta);
    const dtype c = std::cos(0.5 * theta);

    // 0*c + 1*s, where 0 and 1 are logical codewords for Steane code!
    std::vector<dtype> output = {c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0)};
    
    return output;
}

std::vector<dtype> QEC::getLogicalDensityMatrix(const dtype theta) {
    std::vector<dtype> cArrCodeword(FastPauli::CArrLen<NUMQUBITS>);
    const std::vector<dtype> ketL = QEC::getLogicalKet(theta);
    FastPauli::ketToCArr(cArrCodeword, ketL, NUMQUBITS);
    return cArrCodeword;
}
