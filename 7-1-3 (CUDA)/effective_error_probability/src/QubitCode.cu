// Has core circuit and QECC stuff!!!

#include "QubitCode.cuh"

Eigen::Matrix2f QubitCode::getPauli(const char p){
    switch (p){
        case 'I':
            return I_h;
        case 'X':
            return X_h;
        case 'Y':
            return Y_h;
        case 'Z':
            return Z_h;
        default:
            return I_h;
    }
}

void QubitCode::Form_Satbilizers(const std::string &filename){
    Eigen::MatrixXf temp(1, 1);
    std::string stabilizer_line;

    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    while(file >> stabilizer_line){
        temp(0, 0) = {1.0};

        for(auto mat : stabilizer_line){
            temp = qpp::kron(temp, getPauli(mat));
            // temp = qpp::kron(temp, pauli_mat.at(mat));      // Wanted to use unordered_map but apparently it's slower than switch for char/small int cases
        }

        Stabilizers.push_back(temp);
    }
}

// MUST BE CHANGED FOR EACH QECC
void QubitCode::S1Sequential(std::vector<double> &cArr, const double err, const double errTeleport)
{
    // IIIXXXX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 3, err);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::S2Sequential(std::vector<double> &cArr, const double errTeleport)
{
    // IXXIIXX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 1, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::S3Sequential(std::vector<double> &cArr, const double errTeleport)
{
    // XIXIXIX
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCX_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::S4Sequential(std::vector<double> &cArr, const double errTeleport)
{
    // IIIZZZZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 3, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::S5Sequential(std::vector<double> &cArr, const double errTeleport)
{
    // IZZIIZZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 1, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 5, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

void QubitCode::S6Sequential(std::vector<double> &cArr, const double errTeleport)
{
    // ZIZIZIZ
    constexpr int numQubits = NUMQUBITS + 1;
    const int first_qubit = 1;
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 2, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 4, errTeleport);
    NoisyPauli::NCZ_r(cArr, numQubits, 0, first_qubit + 6, errTeleport);
}

// Can probably do this using a map too
void QubitCode::Sn_Sequential(int n, std::vector<double> &cArr, const double err, const double errTeleport)
{
    switch (n)
    {
    case 1:
        QubitCode::S1Sequential(cArr, err, errTeleport);
        break;
    case 2:
        QubitCode::S2Sequential(cArr, errTeleport);
        break;
    case 3:
        QubitCode::S3Sequential(cArr, errTeleport);
        break;
    case 4:
        QubitCode::S4Sequential(cArr, errTeleport);
        break;
    case 5:
        QubitCode::S5Sequential(cArr, errTeleport);
        break;
    case 6:
        QubitCode::S6Sequential(cArr, errTeleport);
        break;
    default:
        break;
    }
}

void QubitCode::nested_loops(int ancilla, int &corr_ind, std::vector<std::vector<double>> &cArrN, std::vector<double> &cArr_Nplus1, std::vector<double> &result, const std::vector<std::vector<double>> &ketH, const std::vector<double> &cArrH1P, const double err, const double errTeleport, const QubitCode::correctionMode correctionMode, const int qblost_or_flag)
{
    FastPauli::dmKron(cArr_Nplus1, cArrH1P, cArrN[ancilla - 1], 1, NUMQUBITS);
    QubitCode::Sn_Sequential(ancilla, cArr_Nplus1, err, errTeleport);
    // Binary
    for (int i = 0; i < 2; i++)
    {
        QubitCode::applyAncillaKet(cArrN[ancilla], 1, NUMQUBITS, ketH[i], cArr_Nplus1);

        if (ancilla == NUM_ANCILLA)
        {
            switch (correctionMode)
            {
            case QubitCode::correctionMode::noLossW1:

                // QubitCode::applyCorrection(corr_ind, cArrN[NUM_ANCILLA]);

                if(corr_ind) QubitCode::applyCorrection(corr_ind - 1, cArrN[NUM_ANCILLA]);

                std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<double>());

                break;
            case QubitCode::correctionMode::oneLoss:

                // QubitCode::applyCorrection_1erasure(corr_ind, cArrN[NUM_ANCILLA], qblost_or_flag);

                if(corr_ind) QubitCode::applyCorrection_1erasure(corr_ind - 1, cArrN[NUM_ANCILLA], qblost_or_flag);

                std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<double>());

                break;
            case QubitCode::correctionMode::noLossW2:

                // QubitCode::applyCorrection_flag(corr_ind, qblost_or_flag, cArrN[NUM_ANCILLA]);

                if(corr_ind) QubitCode::applyCorrection_flag(corr_ind - 1, qblost_or_flag, cArrN[NUM_ANCILLA]);

                std::transform(result.begin(), result.end(), cArrN[NUM_ANCILLA].begin(), result.begin(), std::plus<double>());

                break;

            default:
                break;
            }

            corr_ind++;
        }
        else
        {
            nested_loops(ancilla + 1, corr_ind, cArrN, cArr_Nplus1, result, ketH, cArrH1P, err, errTeleport, correctionMode, qblost_or_flag);
        }
    }
}

// THIS IS BEING USED IN MAIN
std::vector<double> QubitCode::OneRoundSequential(const std::vector<double> &cArrCodeword, const double err, const double errTeleport, const QubitCode::correctionMode correctionMode, const int qblost_or_flag)
{
    const double num = 1 / std::sqrt(2.0);
    // Get ready |+⟩.
    std::vector<double> ketH1P = {num, num};

    // Get ready |+⟩⟨+|.
    std::vector<double> cArrH1P(FastPauli::CArrLen<1>);
    FastPauli::ketToCArr(cArrH1P, ketH1P, 1);

    // Get ready |-⟩.
    std::vector<double> ketH1N = {num, -num};

    std::vector<std::vector<double>> ketH = {ketH1P, ketH1N};

    std::vector<std::vector<double>> cArrN(NUM_ANCILLA + 1, std::vector<double>(FastPauli::CArrLen<NUMQUBITS>));
    std::vector<double> cArr_Nplus1(FastPauli::CArrLen<NUMQUBITS + 1>, 0.0);
    std::vector<double> result(FastPauli::CArrLen<NUMQUBITS>, 0.0);

    // Copy cArrCodeword to cArrN[0]
    cArrN[0] = cArrCodeword;

    int ancilla = 1, corr_ind = 0;

    nested_loops(ancilla, corr_ind, cArrN, cArr_Nplus1, result, ketH, cArrH1P, err, errTeleport, correctionMode, qblost_or_flag);

    return result;
}

void QubitCode::applyCorrection(const uint16_t i, std::vector<double> &cArr)
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

// See QubitCode.cuh for explanation
void QubitCode::applyAncillaKet(std::vector<double> &out_cArr, const int numAncillaQubits, const int numQubits, const std::vector<double> &ancillaKet, const std::vector<double> &in_cArr)
{
    const int nA = 1 << numAncillaQubits;
    const int nQ = 1 << numQubits;
    const int n = 1 << (numQubits + numAncillaQubits);
    int cIndex, cIndex2, cIndex3, cIndex4;
    double sum;
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