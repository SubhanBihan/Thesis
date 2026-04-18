#include "npauli.hpp"
#include "fpauli.hpp"
#include <cstddef>
#include <cstdint>

void NoisyPauli::SingleQubitDepolarise_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit, const dtype err) {
    const uint64_t cArrLen = FastPauli::getCArrLen(numQubits);
    std::vector<dtype> x = cArr;
    std::vector<dtype> y = cArr;
    std::vector<dtype> z = cArr;

    FastPauli::X_r(x, numQubits, targetQubit);
    FastPauli::Y_r(y, numQubits, targetQubit);
    FastPauli::Z_r(z, numQubits, targetQubit);

    for (int i = 0; i < cArrLen; i++) {
        cArr[i] = (1 - err) * cArr[i] + (err / 3.0) * (x[i] + y[i] + z[i]);
    }
}

void NoisyPauli::NCX_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const dtype err) {
    FastPauli::CX_r(cArr, numQubits, controlQubit, targetQubit);

    NoisyPauli::TwoQubitDepolarise_r(cArr, numQubits, controlQubit, targetQubit, err);
}

void NoisyPauli::NCZ_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const dtype err) {
    FastPauli::CZ_r(cArr, numQubits, controlQubit, targetQubit);

    NoisyPauli::TwoQubitDepolarise_r(cArr, numQubits, controlQubit, targetQubit, err);
}

void NoisyPauli::TwoQubitDepolarise_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const dtype err) {
    const uint64_t cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!

    std::vector<dtype> i1x2, i1y2, i1z2, x1x2, x1y2, x1z2, x1i2, y1x2, y1y2, y1z2, y1i2, z1x2, z1y2, z1z2, z1i2;

    i1x2 = cArr;
    i1y2 = cArr;
    i1z2 = cArr;
    x1i2 = cArr;
    y1i2 = cArr;
    z1i2 = cArr;

    FastPauli::X_r(i1x2, numQubits, targetQubit);
    x1x2 = i1x2;
    y1x2 = i1x2;
    z1x2 = i1x2;
    FastPauli::X_r(x1x2, numQubits, controlQubit);
    FastPauli::Y_r(y1x2, numQubits, controlQubit);
    FastPauli::Z_r(z1x2, numQubits, controlQubit);

    FastPauli::Y_r(i1y2, numQubits, targetQubit);
    x1y2 = i1y2;
    y1y2 = i1y2;
    z1y2 = i1y2;
    FastPauli::X_r(x1y2, numQubits, controlQubit);
    FastPauli::Y_r(y1y2, numQubits, controlQubit);
    FastPauli::Z_r(z1y2, numQubits, controlQubit);

    FastPauli::Z_r(i1z2, numQubits, targetQubit);
    x1z2 = i1z2;
    y1z2 = i1z2;
    z1z2 = i1z2;
    FastPauli::X_r(x1z2, numQubits, controlQubit);
    FastPauli::Y_r(y1z2, numQubits, controlQubit);
    FastPauli::Z_r(z1z2, numQubits, controlQubit);

    for (size_t i = 0; i < cArrLen; i++) {
        cArr[i] = (1 - err) * cArr[i] + (err / 15.0) * (i1x2[i] +
                                                        i1y2[i] +
                                                        i1z2[i] +
                                                        x1x2[i] +
                                                        x1y2[i] +
                                                        x1z2[i] +
                                                        x1i2[i] +
                                                        y1x2[i] +
                                                        y1y2[i] +
                                                        y1z2[i] +
                                                        y1i2[i] +
                                                        z1x2[i] +
                                                        z1y2[i] +
                                                        z1z2[i] +
                                                        z1i2[i]);
    }
}