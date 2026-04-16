// Need to change encoded logical qubit (output of getLogicalKet)
// Changes have been made

#include "FiveQubitCode.hpp"

std::vector<dtype> FiveQubitCode::getLogicalKet(const dtype theta) {
    const dtype s = static_cast<dtype>(std::sin(0.5 * theta));
    const dtype c = static_cast<dtype>(std::cos(0.5 * theta));

    // 0*c + 1*s, where 0 and 1 are logical codewords for Steane code!
    return {0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, -0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * c, 0.0, 0.0, -0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * s, 0.0, 0.0, 0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * c, 0.0, 0.0, 0.0, 0.0, -0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25f * s, -0.25f * c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25f * c};
}

std::vector<dtype> FiveQubitCode::getLogicalDensityMatrix(const dtype theta) {
    std::vector<dtype> cArrCodeword(FastPauli::CArrLen<NUMQUBITS>);
    const std::vector<dtype> ketL = FiveQubitCode::getLogicalKet(theta);
    FastPauli::ketToCArr(cArrCodeword, ketL, NUMQUBITS);
    return cArrCodeword;
}
