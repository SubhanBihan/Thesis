// Has core circuit and QECC stuff!!!

#include "QEC.hpp"
// #include <stack>

std::vector<dtype> ketH1P = {num, num};

// Get ready |+⟩⟨+|.
std::vector<dtype> cArrH1P(FastPauli::CArrLen<1>);

// Get ready |-⟩.
std::vector<dtype> ketH1N = {num, -num};

std::vector<std::vector<dtype>> ketH = {ketH1P, ketH1N};

// MUST BE CHANGED FOR EACH QECC
void QEC::S1Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // IIIXXXX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 3, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S2Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // IXXIIXX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 1, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S3Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // XIXIXIX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S4Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // IIIZZZZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 3, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, err);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S5Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // IZZIIZZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 1, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, err);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QEC::S6Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    // ZIZIZIZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, err);
}

void QEC::Sn_Sequential(int n, std::vector<dtype> &cArr, const dtype err, const dtype errTeleport)
{
    switch (n)
    {
        case 1:
            QEC::S1Sequential(cArr, err, errTeleport);
            break;
        case 2:
            QEC::S2Sequential(cArr, err, errTeleport);
            break;
        case 3:
            QEC::S3Sequential(cArr, err, errTeleport);
            break;
        case 4:
            QEC::S4Sequential(cArr, err, errTeleport);
            break;
        case 5:
            QEC::S5Sequential(cArr, err, errTeleport);
            break;
        case 6:
            QEC::S6Sequential(cArr, err, errTeleport);
            break;
        default:
            break;
    }
}

void QEC::ORS_NestedLoops(int ancilla, int &corr_ind, std::vector<std::vector<dtype>> &cArrN, std::vector<dtype> &result, const dtype err, const dtype errTeleport, const QEC::correctionMode correctionMode, const int qblost_or_flag)
{
    std::vector<dtype> cArr_Nplus1(FastPauli::CArrLen<NUMQUBITS + 1>);

    FastPauli::dmKron(cArr_Nplus1, cArrH1P, cArrN[ancilla - 1], 1, NUMQUBITS);
    QEC::Sn_Sequential(ancilla, cArr_Nplus1, err, errTeleport);

    for (int i = 0; i < 2; i++)
    {
        QEC::applyAncillaKet(cArrN[ancilla], 1, NUMQUBITS, ketH[i], cArr_Nplus1);

        if (ancilla == NUM_ANCILLA)
        {
            if(corr_ind)
                switch (correctionMode)
                {
                    case QEC::correctionMode::noLossW1:
                        QEC::applyCorrection(corr_ind - 1, cArrN[NUM_ANCILLA]);
                        break;
                    case QEC::correctionMode::oneLoss:
                        QEC::applyCorrection_1erasure(corr_ind - 1, cArrN[NUM_ANCILLA], qblost_or_flag);
                        break;
                    case QEC::correctionMode::noLossW2:
                        QEC::applyCorrection_flag(corr_ind - 1, qblost_or_flag, cArrN[NUM_ANCILLA]);
                        break;
                }

            std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<dtype>());

            corr_ind++;
        }
        else
        {
            ORS_NestedLoops(ancilla + 1, corr_ind, cArrN, result, err, errTeleport, correctionMode, qblost_or_flag);
        }
    }
}

// THIS IS BEING USED IN MAIN
std::vector<dtype> QEC::OneRoundSequential(const std::vector<dtype> &cArrCodeword, const dtype err, const dtype errTeleport, const QEC::correctionMode correctionMode, const int qblost_or_flag)
{
    std::vector<std::vector<dtype>> cArrN(NUM_ANCILLA + 1, std::vector<dtype>(FastPauli::CArrLen<NUMQUBITS>));
    std::vector<dtype> result(FastPauli::CArrLen<NUMQUBITS>, 0.0);

    // Copy cArrCodeword to cArrN[0]
    cArrN[0] = cArrCodeword;

    int ancilla = 1, corr_ind = 0;

    ORS_NestedLoops(ancilla, corr_ind, cArrN, result, err, errTeleport, correctionMode, qblost_or_flag);

    return result;
}

// // Iterative DFS approach (Equivalent)
// std::vector<dtype> QEC::OneRoundSequential(const std::vector<dtype> &cArrCodeword, const dtype err, const dtype errTeleport, const QEC::correctionMode correctionMode, const int qblost_or_flag)
// {
//     std::vector<dtype> result(FastPauli::CArrLen<NUMQUBITS>, 0.0);
//     // Stack stores the ancilla level and the state vector
//     std::stack<std::pair<int, std::vector<dtype>>> stk;
//     stk.push({1, cArrCodeword});
    
//     int corr_ind = 0;

//     while (!stk.empty())
//     {
//         auto [ancilla, current_cArr] = std::move(stk.top());
//         stk.pop();

//         std::vector<dtype> cArr_Nplus1(FastPauli::CArrLen<NUMQUBITS + 1>);
//         FastPauli::dmKron(cArr_Nplus1, cArrH1P, current_cArr, 1, NUMQUBITS);
//         QEC::Sn_Sequential(ancilla, cArr_Nplus1, err, errTeleport);

//         if (ancilla == NUM_ANCILLA)
//         {
//             // 2. LEAF LOGIC: Must process i = 0 then i = 1 to match corr_ind order
//             for (int i = 0; i < 2; i++)
//             {
//                 std::vector<dtype> cArr_leaf(FastPauli::CArrLen<NUMQUBITS>);
//                 QEC::applyAncillaKet(cArr_leaf, 1, NUMQUBITS, ketH[i], cArr_Nplus1);

//                 if (corr_ind > 0) { // Equivalent to if(corr_ind)
//                     switch (correctionMode)
//                     {
//                         case QEC::correctionMode::noLossW1:
//                             QEC::applyCorrection(corr_ind - 1, cArr_leaf); break;
//                         case QEC::correctionMode::oneLoss:
//                             QEC::applyCorrection_1erasure(corr_ind - 1, cArr_leaf, qblost_or_flag); break;
//                         case QEC::correctionMode::noLossW2:
//                             QEC::applyCorrection_flag(corr_ind - 1, qblost_or_flag, cArr_leaf); break;
//                     }
//                 }
                
//                 std::transform(result.begin(), result.end(), cArr_leaf.begin(), result.begin(), std::plus<dtype>());
//                 corr_ind++;
//             }
//         }
//         else
//         {
//             // 3. INTERNAL NODES: Push i = 1 then i = 0 so that i = 0 is popped first
//             for (int i = 1; i >= 0; i--)
//             {
//                 std::vector<dtype> next_cArr(FastPauli::CArrLen<NUMQUBITS>);
//                 QEC::applyAncillaKet(next_cArr, 1, NUMQUBITS, ketH[i], cArr_Nplus1);
//                 stk.push({ancilla + 1, std::move(next_cArr)});
//             }
//         }
//     }
//     return result;
// }

void QEC::applyCorrection(const uint16_t i, std::vector<dtype> &cArr)
{
    for (int j = 0; j < NUMQUBITS; j++)
    {
        switch (corrections[i][j])
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

// See QEC.hpp for explanation
void QEC::applyAncillaKet(std::vector<dtype> &out_cArr, const int numAncillaQubits, const int numQubits, const std::vector<dtype> &ancillaKet, const std::vector<dtype> &in_cArr)
{
    const int nA = 1 << numAncillaQubits;
    const int nQ = 1 << numQubits;
    const int n = 1 << (numQubits + numAncillaQubits);

    for (int p = 0; p < nQ; p++)
    {
        for (int q = p; q < nQ; q++)
        {
            dtype sum = 0;
            for (int i = 0; i < nA; i++)
            {
                sum += in_cArr[FastPauli::getRow(i * nQ + p, i * nQ + q, n)] * ancillaKet[i] * ancillaKet[i];
                for (int j = i + 1; j < nA; j++)
                {
                    sum += (in_cArr[FastPauli::getRow(i * nQ + p, j * nQ + q, n)] + in_cArr[FastPauli::getRow(i * nQ + q, j * nQ + p, n)]) * ancillaKet[i] * ancillaKet[j];
                }
            }
            out_cArr[FastPauli::getRow(p, q, nQ)] = sum;
        }
    }
}