// Has core circuit and QECC stuff!!!

#include "FiveQubitCode.hpp"

std::vector<dtype> ketH1P(2);

// Get ready |+⟩⟨+|.
std::vector<dtype> cArrH1P(FastPauli::CArrLen<1>);

// Get ready |-⟩.
std::vector<dtype> ketH1N(2);

std::vector<std::vector<dtype>> ketH(2);

// MUST BE CHANGED FOR EACH QECC
void FiveQubitCode::S1Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // IIIXXXX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 3, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void FiveQubitCode::S2Sequential(std::vector<dtype> &cArr, const dtype errTeleport)
{
    // IXXIIXX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 1, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void FiveQubitCode::S3Sequential(std::vector<dtype> &cArr, const dtype errTeleport)
{
    // XIXIXIX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void FiveQubitCode::S4Sequential(std::vector<dtype> &cArr, const dtype errTeleport)
{
    // IIIZZZZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 3, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void FiveQubitCode::S5Sequential(std::vector<dtype> &cArr, const dtype errTeleport)
{
    // IZZIIZZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 1, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void FiveQubitCode::S6Sequential(std::vector<dtype> &cArr, const dtype errTeleport)
{
    // ZIZIZIZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void FiveQubitCode::Sn_Sequential(int n, std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    switch (n)
    {
    case 1:
        FiveQubitCode::S1Sequential(cArr, err, errTeleport);
        break;
    case 2:
        FiveQubitCode::S2Sequential(cArr, errTeleport);
        break;
    case 3:
        FiveQubitCode::S3Sequential(cArr, errTeleport);
        break;
    case 4:
        FiveQubitCode::S4Sequential(cArr, errTeleport);
        break;
    case 5:
        FiveQubitCode::S5Sequential(cArr, errTeleport);
        break;
    case 6:
        FiveQubitCode::S6Sequential(cArr, errTeleport);
        break;
    default:
        break;
    }
}

void FiveQubitCode::ORS_NestedLoops(int ancilla, int &corr_ind, std::vector<std::vector<dtype>> &cArrN, std::vector<dtype> &result, const dtype err, const dtype errTeleport, const FiveQubitCode::correctionMode correctionMode, const int qblost_or_flag)
{
    std::vector<dtype> cArr_Nplus1(FastPauli::CArrLen<NUMQUBITS + 1>);

    FastPauli::dmKron(cArr_Nplus1, cArrH1P, cArrN[ancilla - 1], 1, NUMQUBITS);
    FiveQubitCode::Sn_Sequential(ancilla, cArr_Nplus1, err, errTeleport);

    for (int i = 0; i < 2; i++)
    {
        FiveQubitCode::applyAncillaKet(cArrN[ancilla], 1, NUMQUBITS, ketH[i], cArr_Nplus1);

        if (ancilla == NUM_ANCILLA)
        {
            switch (correctionMode)
            {
            case FiveQubitCode::correctionMode::noLossW1:
                if(corr_ind) FiveQubitCode::applyCorrection(corr_ind - 1, cArrN[NUM_ANCILLA]);

                std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<dtype>());

                break;

            case FiveQubitCode::correctionMode::oneLoss:
                if(corr_ind) FiveQubitCode::applyCorrection_1erasure(corr_ind - 1, cArrN[NUM_ANCILLA], qblost_or_flag);

                std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<dtype>());

                break;

            case FiveQubitCode::correctionMode::noLossW2:
                if(corr_ind) FiveQubitCode::applyCorrection_flag(corr_ind - 1, qblost_or_flag, cArrN[NUM_ANCILLA]);

                std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<dtype>());

                break;

            default:
                break;
            }

            corr_ind++;
        }
        else
        {
            ORS_NestedLoops(ancilla + 1, corr_ind, cArrN, result, err, errTeleport, correctionMode, qblost_or_flag);
        }
    }
}

// THIS IS BEING USED IN MAIN
std::vector<dtype> FiveQubitCode::OneRoundSequential(const std::vector<dtype> &cArrCodeword, const dtype err, const dtype errTeleport, const FiveQubitCode::correctionMode correctionMode, const int qblost_or_flag)
{
    std::vector<std::vector<dtype>> cArrN(NUM_ANCILLA + 1, std::vector<dtype>(FastPauli::CArrLen<NUMQUBITS>));

    std::vector<dtype> result(FastPauli::CArrLen<NUMQUBITS>, 0.0);

    // Copy cArrCodeword to cArrN[0]
    cArrN[0] = cArrCodeword;

    int ancilla = 1, corr_ind = 0;

    ORS_NestedLoops(ancilla, corr_ind, cArrN, result, err, errTeleport, correctionMode, qblost_or_flag);

    return result;
}

void FiveQubitCode::applyCorrection(const uint16_t i, std::vector<dtype> &cArr)
{
    // if (!i)
    //     return;
        
    // std::string correction_i = corrections[i - 1];
    std::string correction_i = corrections[i];

    for (int j = 0; j < NUMQUBITS; j++)
    {
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

// See FiveQubitCode.hpp for explanation
void FiveQubitCode::applyAncillaKet(std::vector<dtype> &out_cArr, const int numAncillaQubits, const int numQubits, const std::vector<dtype> &ancillaKet, const std::vector<dtype> &in_cArr)
{
    const int nA = 1 << numAncillaQubits;
    const int nQ = 1 << numQubits;
    const int n = 1 << (numQubits + numAncillaQubits);
    int cIndex, cIndex2, cIndex3, cIndex4;
    dtype sum;
    for (int p = 0; p < nQ; p++)
    {
        for (int q = p; q < nQ; q++)
        {
            cIndex = FastPauli::getRow(p, q, nQ);
            sum = 0;
            for (int i = 0; i < nA; i++)
            {
                cIndex2 = FastPauli::getRow(i * nQ + p, i * nQ + q, n);
                sum += in_cArr[cIndex2] * ancillaKet[i] * ancillaKet[i];
                for (int j = i + 1; j < nA; j++)
                {
                    cIndex3 = FastPauli::getRow(i * nQ + p, j * nQ + q, n);
                    cIndex4 = FastPauli::getRow(i * nQ + q, j * nQ + p, n);
                    sum += (in_cArr[cIndex3] + in_cArr[cIndex4]) * ancillaKet[i] * ancillaKet[j];
                }
            }
            out_cArr[cIndex] = sum;
        }
    }
}