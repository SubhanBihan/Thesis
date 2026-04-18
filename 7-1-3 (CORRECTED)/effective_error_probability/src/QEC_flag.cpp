// Has core circuit and QECC stuff!!! Must edit extensively!

#include "QEC.hpp"

// MUST BE CHANGED FOR EACH QECC
// For these Sn_flag functions, ancilla and flag are at 0 and 1 index resp. 
// All but 1 of the ancilla-data gates will have errTeleport
void QEC::S1_flag(std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // IIIXXXX
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 3, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S2_flag(std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // IXXIIXX
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 1, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S3_flag(std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // XIXIXIX
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S4_flag(std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // IIIZZZZ
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 3, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, err);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S5_flag(std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // IZZIIZZ
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 1, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S6_flag(std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // ZIZIZIZ
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, err);
}

void QEC::Sn_flag(int n, std::vector<dtype>& cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr){
    switch (n)
    {
        case 1:
            QEC::S1_flag(cArr, err, errTeleport, ancillaFlagErr);
            break;
        case 2:
            QEC::S2_flag(cArr, err, errTeleport, ancillaFlagErr);
            break;
        case 3:
            QEC::S3_flag(cArr, err, errTeleport, ancillaFlagErr);
            break;
        case 4:
            QEC::S4_flag(cArr, err, errTeleport, ancillaFlagErr);
            break;
        case 5:
            QEC::S5_flag(cArr, err, errTeleport, ancillaFlagErr);
            break;
        case 6:
            QEC::S6_flag(cArr, err, errTeleport, ancillaFlagErr);
            break;
        default:
            break;
    }
}

// Seems to check for both untriggered ancilla and flag qubit (1 each): will only output data qubits if so
std::vector<dtype> QEC::getcArrW0(const std::vector<dtype>& cArr) {
    std::vector<dtype> Ket_AF_00 = QEC::getAncillaFlagKet<false, false>();
    std::vector<dtype> result = QEC::applyAncillaKet<2, NUMQUBITS>(Ket_AF_00, cArr);
    return result;
}

// Check for ancilla triggered, flag not triggered
std::vector<dtype> QEC::getcArrW1(const std::vector<dtype>& cArr) {
    std::vector<dtype> Ket_AF_10 = QEC::getAncillaFlagKet<true, false>();
    std::vector<dtype> result = QEC::applyAncillaKet<2, NUMQUBITS>(Ket_AF_10, cArr);
    return result;
}

// Checks for flag triggered, ancilla both triggered and not triggered
// The resulting data qubits from these two are added - well, only one of them will be non-zero
// Essentially we get the data qubits ONLY IF flag is triggered
std::vector<dtype> QEC::getcArrW2(const std::vector<dtype>& cArr) {
    std::vector<dtype> Ket_AF_01 = QEC::getAncillaFlagKet<false, true>();
    std::vector<dtype> Ket_AF_11 = QEC::getAncillaFlagKet<true, true>();
    std::vector<dtype> tmp = QEC::applyAncillaKet<2, NUMQUBITS>(Ket_AF_11, cArr);
    std::vector<dtype> result = QEC::applyAncillaKet<2, NUMQUBITS>(Ket_AF_01, cArr);
    std::transform(result.begin(), result.end(), tmp.begin(), result.begin(), std::plus<dtype>());
    return result;
}

std::vector<dtype> QEC::OneRoundSequential_flag(const std::vector<dtype> &logicalcArr, const std::vector<dtype> &cArrCodeword, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr) {
    // std::cout << "START OneRoundSequential_flag" << std::endl;
    // Initialise the ancilla-flag pair density matrix. Only the untriggered one is needed.
    std::vector<dtype> AF_00 = QEC::getAncillaFlagDensityMatrix<false, false>();
    // AuxMath::print1DVec(AF_00);

    // Initialise temporary kets for performing intermediary calculations.
    std::vector<dtype> tmp(FastPauli::CArrLen<NUMQUBITS + 2>, 0.0);
    // std::cout << "After tmp array alloc" << std::endl;
    std::vector<std::vector<dtype>> cArrN(3, std::vector<dtype>(FastPauli::CArrLen<NUMQUBITS>));
    cArrN[0] = cArrCodeword;
    // AuxMath::print1DVec(cArrN[0]);
    // std::cout << "After BIG array alloc, before loop" << std::endl;
    
    std::vector<dtype> result(FastPauli::CArrLen<NUMQUBITS>, 0.0);

    // This works according to Definition 3.4.2 in KJW's thesis
    // ALSO Aggregate all corrected matrices. Use cArrW2 to store the results to avoid declaring new std::vector.
    // ONLY 1 result is acutally non-zero
    for(int ancilla = 1; ancilla <= NUM_ANCILLA; ancilla++){
        // std::cout << "ancilla = " << ancilla << std::endl;
        tmp = FastPauli::dmKron<2, NUMQUBITS>(AF_00, cArrN[0]);
        // std::cout << "1" << std::endl;
        QEC::Sn_flag(ancilla, tmp, err, errTeleport, ancillaFlagErr);
        // std::cout << "2" << std::endl;
        cArrN[0] = QEC::getcArrW0(tmp);
        // std::cout << "3" << std::endl;
        cArrN[1] = QEC::getcArrW1(tmp);
        // std::cout << "4" << std::endl;
        cArrN[2] = QEC::getcArrW2(tmp);
        // std::cout << "5" << std::endl;

        cArrN[1] = QEC::OneRoundSequential(cArrN[1], err, errTeleport, QEC::correctionMode::noLossW1, 0);
        // std::cout << "6" << std::endl;
        cArrN[2] = QEC::OneRoundSequential(cArrN[2], err, errTeleport, QEC::correctionMode::noLossW2, ancilla - 1);
        // std::cout << "7" << std::endl;

        // std::transform(result.begin(), result.end(), cArrN[1].begin(), result.begin(), std::plus<dtype>());
        // std::transform(result.begin(), result.end(), cArrN[2].begin(), result.begin(), std::plus<dtype>());
        for (size_t i = 0; i < result.size(); ++i)
            result[i] += cArrN[1][i] + cArrN[2][i];
        // std::cout << "8" << std::endl;
    }

    // std::cout << "Loop ended" << std::endl;

    std::transform(result.begin(), result.end(), cArrN[0].begin(), result.begin(), std::plus<dtype>());
    
    // std::cout << "After transform" << std::endl;

    return result;
}

void QEC::applyCorrection_flag(const uint16_t i, const uint16_t flag, std::vector<dtype>& cArr) {
    for(int j = 0; j < NUMQUBITS; j++){
        switch (corrections_flag[flag * (DIM_ANCILLA - 1) + i][j])
        {
            case 'I':
                continue;
            case 'X':
                FastPauli::X_r(cArr, NUMQUBITS, j);
                break;
            case 'Y':
                FastPauli::Y_r(cArr, NUMQUBITS, j);
                break;
            case 'Z':
                FastPauli::Z_r(cArr, NUMQUBITS, j);
                break;
            default:
                break;
        }
    }
}