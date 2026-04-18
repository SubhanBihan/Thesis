#include "FiveQubitCode.hpp"

void FiveQubitCode::applyCorrection_1erasure(const uint16_t i, std::vector<dtype>& cArr, const int k2)
{   
    // k2 is qubitlostindex
    for(int j = 0; j < NUMQUBITS; j++){
        switch (corrections_erasure[k2 * (DIM_ANCILLA - 1) + i][j])
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