// Has core circuit and QECC stuff!!! Must edit extensively!

#include "QubitCode.cuh"

// MUST BE CHANGED FOR EACH QECC
// For these Sn_flag functions, ancilla and flag are at 0 and 1 index resp. 
// All but 1 of the ancilla-data gates will have errTeleport
void QubitCode::S1_flag(std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr) {
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

void QubitCode::S2_flag(std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr) {
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

void QubitCode::S3_flag(std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr) {
    // XIXIXIX
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::S4_flag(std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr) {
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

void QubitCode::S5_flag(std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr) {
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

void QubitCode::S6_flag(std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr) {
    // ZIZIZIZ
    const int first_qubit = 2;
    constexpr int numQubits = NUMQUBITS + 2;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, err);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, 1, ancillaFlagErr);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::Sn_flag(int n, std::vector<double>& cArr, const double err, const double errTeleport, const double ancillaFlagErr){
    switch (n)
    {
    case 1:
        QubitCode::S1_flag(cArr, err, errTeleport, ancillaFlagErr);
        break;
    case 2:
        QubitCode::S2_flag(cArr, err, errTeleport, ancillaFlagErr);
        break;
    case 3:
        QubitCode::S3_flag(cArr, err, errTeleport, ancillaFlagErr);
        break;
    case 4:
        QubitCode::S4_flag(cArr, err, errTeleport, ancillaFlagErr);
        break;
    case 5:
        QubitCode::S5_flag(cArr, err, errTeleport, ancillaFlagErr);
        break;
    case 6:
        QubitCode::S6_flag(cArr, err, errTeleport, ancillaFlagErr);
        break;
    default:
        break;
    }
}

// Seems to check for both untriggered ancilla and flag qubit (1 each): will only output data qubits if so
std::vector<double> QubitCode::getcArrW0(const std::vector<double>& cArr) {
    std::vector<double> Ket_AF_00 = QubitCode::getAncillaFlagKet<false, false>();
    std::vector<double> result = QubitCode::applyAncillaKet<2, NUMQUBITS>(Ket_AF_00, cArr);
    return result;
}

// Check for ancilla triggered, flag not triggered
std::vector<double> QubitCode::getcArrW1(const std::vector<double>& cArr) {
    std::vector<double> Ket_AF_10 = QubitCode::getAncillaFlagKet<true, false>();
    std::vector<double> result = QubitCode::applyAncillaKet<2, NUMQUBITS>(Ket_AF_10, cArr);
    return result;
}

// Checks for flag triggered, ancilla both triggered and not triggered
// The resulting data qubits from these two are added - well, only one of them will be non-zero
// Essentially we get the data qubits ONLY IF flag is triggered
std::vector<double> QubitCode::getcArrW2(const std::vector<double>& cArr) {
    std::vector<double> Ket_AF_01 = QubitCode::getAncillaFlagKet<false, true>();
    std::vector<double> Ket_AF_11 = QubitCode::getAncillaFlagKet<true, true>();
    std::vector<double> tmp = QubitCode::applyAncillaKet<2, NUMQUBITS>(Ket_AF_11, cArr);
    std::vector<double> result = QubitCode::applyAncillaKet<2, NUMQUBITS>(Ket_AF_01, cArr);
    std::transform(result.begin(), result.end(), tmp.begin(), result.begin(), std::plus<double>());
    return result;
}

std::vector<double> QubitCode::OneRoundSequential_flag(const std::vector<double> &logicalcArr, const std::vector<double> &cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr) {
    // std::cout << "START OneRoundSequential_flag" << std::endl;
    // Initialise the ancilla-flag pair density matrix. Only the untriggered one is needed.
    std::vector<double> AF_00 = QubitCode::getAncillaFlagDensityMatrix<false, false>();
    // AuxMath::print1DVec(AF_00);

    // Initialise temporary kets for performing intermediary calculations.
    std::vector<double> tmp(FastPauli::CArrLen<NUMQUBITS + 2>, 0.0);
    // std::cout << "After tmp array alloc" << std::endl;
    std::vector<std::vector<double>> cArrN(3, std::vector<double>(FastPauli::CArrLen<NUMQUBITS>));
    cArrN[0] = cArrCodeword;
    // AuxMath::print1DVec(cArrN[0]);
    // std::cout << "After BIG array alloc, before loop" << std::endl;
    
    std::vector<double> result(FastPauli::CArrLen<NUMQUBITS>, 0.0);

    // This works according to Definition 3.4.2 in KJW's thesis
    // ALSO Aggregate all corrected matrices.
    // ONLY 1 result is acutally non-zero
    for(int ancilla = 1; ancilla <= NUM_ANCILLA; ancilla++){
        // std::cout << "ancilla = " << ancilla << std::endl;
        tmp = FastPauli::dmKron<2, NUMQUBITS>(AF_00, cArrN[0]);
        // std::cout << "1" << std::endl;
        QubitCode::Sn_flag(ancilla, tmp, err, errTeleport, ancillaFlagErr);
        // std::cout << "2" << std::endl;
        cArrN[0] = QubitCode::getcArrW0(tmp);
        // std::cout << "3" << std::endl;
        cArrN[1] = QubitCode::getcArrW1(tmp);
        // std::cout << "4" << std::endl;
        cArrN[2] = QubitCode::getcArrW2(tmp);
        // std::cout << "5" << std::endl;

        cArrN[1] = QubitCode::OneRoundSequential(cArrN[1], err, errTeleport, QubitCode::correctionMode::noLossW1, 0);
        // std::cout << "6" << std::endl;
        cArrN[2] = QubitCode::OneRoundSequential(cArrN[2], err, errTeleport, QubitCode::correctionMode::noLossW2, ancilla - 1);
        // std::cout << "7" << std::endl;

        std::transform(result.begin(), result.end(), cArrN[1].begin(), result.begin(), std::plus<double>());
        std::transform(result.begin(), result.end(), cArrN[2].begin(), result.begin(), std::plus<double>());
        // std::cout << "8" << std::endl;
    }

    // std::cout << "Loop ended" << std::endl;

    std::transform(result.begin(), result.end(), cArrN[0].begin(), result.begin(), std::plus<double>());
    
    // std::cout << "After transform" << std::endl;

    return result;
}

void QubitCode::applyCorrection_flag(const uint16_t i, const uint16_t flag, std::vector<double>& cArr) {
    // if(!i) return;
    // std::string correction_i = corrections_flag[flag * (DIM_ANCILLA - 1) + i - 1];

    std::string correction_i = corrections_flag[flag * (DIM_ANCILLA - 1) + i];
    for(int j = 0; j < NUMQUBITS; j++){
        switch (correction_i[j])
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