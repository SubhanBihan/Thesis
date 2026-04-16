// Need to change encoded logical qubit (output of getLogicalKet) according to the QECC

// #include <cmath>

#include "QubitCode.cuh"

std::vector<datatype> QubitCode::getLogicalKet(const double theta) {
    const datatype s = std::sin(0.5 * theta);
    const datatype c = std::cos(0.5 * theta);

    auto fsqrt = [](double x) -> float { return static_cast<datatype>(std::sqrt(x)); };

    // 0*c + 1*s, where 0 and 1 are logical codewords for Stabilizer code!
    std::vector<datatype> ket = {c / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / fsqrt(8.0), 0, 0, 0, 0, 0, 0, s / fsqrt(8.0), 0, 0, s / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / fsqrt(8.0), 0, 0, 0, 0, s / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / fsqrt(8.0), 0, 0, 0, 0, 0, 0, s / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, s / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / fsqrt(8.0), 0, 0, 0, 0, c / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / fsqrt(8.0), 0, 0, c / fsqrt(8.0), 0, 0, 0, 0, 0, 0, s / fsqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / fsqrt(8.0)};

    return ket;
}

// CHECK THIS!!!
std::vector<datatype> QubitCode::getLogicalDensityMatrix(const double theta) {
    auto ket = QubitCode::getLogicalKet(theta);
    std::vector<datatype> densityMatrix(ket.size() * ket.size());

    unsigned index = 0;
    // Compute the density matrix from the ket
    for (auto& x : ket) {
        for (auto& y : ket) {
            densityMatrix[index++] = x * y;
        }
    }

    return densityMatrix;
}
