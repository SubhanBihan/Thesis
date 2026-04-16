#include "QEC.hpp"

void QEC::applyCorrection_1erasure(const uint16_t i, std::vector<dtype>& cArr, const int k2)
{   
    // k2 is qubitlostindex
    // ALSO CHECK THE INDEX ISSUE!!! - CHECKED, but need to test properly
    // if(!i) return;
    // std::string correction_i = corrections_erasure[k2 * (DIM_ANCILLA - 1) + i - 1];

    std::string correction_i = corrections_erasure[k2 * (DIM_ANCILLA - 1) + i];
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