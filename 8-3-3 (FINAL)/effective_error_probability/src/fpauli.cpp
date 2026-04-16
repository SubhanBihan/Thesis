#include "fpauli.hpp"

// #define NUMTHREADS 8  // Number of CPU threads

const uint64_t FastPauli::getCArrLen(const int numQubits) {
    return (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));
}

void FastPauli::dmKron(std::vector<dtype>& result, const std::vector<dtype>& CArr1, const std::vector<dtype>& CArr2, int numQubits1, int numQubits2)
{
    // MIGHT WANT TO RESIZE BEFOREHAND
    const int ntot = 1 << (numQubits1 + numQubits2);
    const int n1 = 1 << numQubits1;
    const int n2 = 1 << numQubits2;
    int i, j;
    int tmp, tmp0, tmp1;
    int mi, mj;
    int newInd, oldInd1, oldInd2;
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArrNew, cArrOld1, cArrOld2) private(i, j, mi, mj, tmp, tmp0, tmp1, newInd, oldInd1, oldInd2) schedule(static, std::max(1, ntot / NUMTHREADS))
    // #endif
    for (i = 0; i < ntot; i++) {
        for (j = i; j < ntot; j++) {
            mi = FastPauli::fastMod(i, numQubits2);
            mj = FastPauli::fastMod(j, numQubits2);
            tmp = FastPauli::fastMod(i, numQubits2) - FastPauli::fastMod(j, numQubits2);
            tmp0 = tmp < 0 ? 0 : tmp;
            tmp = FastPauli::fastMod(i - mj, numQubits2);
            tmp1 = tmp > (n2 - mj - 1) ? 0 : tmp;
            newInd = FastPauli::getRow(i, j, ntot);
            oldInd1 = FastPauli::getRow(i / n2, j / n2, n1);
            oldInd2 = FastPauli::getRow(mi - tmp1, mj + tmp0, n2);
            result[newInd] = CArr1[oldInd1] * CArr2[oldInd2];
        }
    }
}

void FastPauli::ketKron(std::vector<dtype>& result, const std::vector<dtype>& ket1, const std::vector<dtype>& ket2, const int numQubits1, const int numQubits2) {
    const int ntot = 1 << (numQubits1 + numQubits2);
    // const int n1 = 1 << numQubits1; // Not used.
    const int n2 = 1 << numQubits2;
    for (int i = 0; i < ntot; i++) {
        result[i] = ket1[i / n2] * ket2[FastPauli::fastMod(i, numQubits2)];
    }
}

void FastPauli::ketToCArr(std::vector<dtype>& cArr, const std::vector<dtype>& ket, const int numQubits) {
    // Computes the ket-bra product matrix - this is symmetric, so only upper-triangular portion is condensed to cArr
    // CArrLen then CArr is thus set according to this rule
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            cArr[FastPauli::getRow(i, j, n)] = ket[i] * ket[j];
        }
    }
}

void FastPauli::apply_permutation(std::vector<dtype>& v, int *indices, const int vSize) {  // Reference: https://stackoverflow.com/questions/67751784/how-to-do-in-place-sorting-a-list-according-to-a-given-index-in-c.
    for (int i = 0; i < vSize; i++) {
        auto current = i;
        while (i != indices[current]) {
            auto next = indices[current];
            std::swap(v[current], v[next]);
            indices[current] = current;
            current = next;
        }
        indices[current] = current;
    }
}

int FastPauli::getRow(int i, int j, const int n) {
    // This function maps a pair of indices i and j (representing rows and columns in an upper triangular matrix) 
    // into a row-wise 1D index for the array cArr[]. It maps the upper triangular part of the matrix into a single array.
    // GO TO LINK to learn more
    // Note: 1 ≤ i, j ≤ n and i ≤ j.
    // Note: n here refers to the dimension of the square matrix, not the number of qubits.
    i++;
    j++;
    assert(1 <= i && j <= n && i <= j);
    return j + n * (i - 1) - ((i * (i - 1)) >> 1) - 1;  // Reference: https://www.geeksforgeeks.org/convert-given-upper-triangular-matrix-to-1d-array/
}

// int getCol(int i, int j){
//     i = i + 1;
//     j = j + 1;
//     // std::cout << "i: " << i << ", j: " << j << "\n";
//     assert(1 <= i && i <= j);
//     // Note: 1 ≤ i, j ≤ n and i ≤ j.
//     // Note: n here refers to the dimension of the square matrix, not the number of qubits.
//     return i + ((j*(j-1)) >> 1) - 1; // Reference: https://www.geeksforgeeks.org/convert-given-upper-triangular-matrix-to-1d-array/
// }

int FastPauli::fastMod(const int i, const int t) {
    return (i) & ((1 << t) - 1);
}

int FastPauli::congMod(const int i, const int t) {
    // Numbers congruent to mod 2^t;
    // In Matlab, it would be '2*i - 2^t + mod(-i,2^t)'. Note that this uses index starting from 1 as opposed to the C/C++ function here which uses index starting from 0.
    // congMod(i,t) gives the numbers that are congruent to {0,1,...,2^t-1} mod 2^(t+1).
    const int pow2t = 1 << t;
    return ((i + 1) << 1) - pow2t + ((-i - 1) & (pow2t - 1)) - 1;  // If b is a power of 2, it is possible to compute n % b as n & (b-1). Reference: https://stackoverflow.com/questions/7594508/modulo-operator-with-negative-values.
}

void FastPauli::ketX_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit) {
    int i;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int tmp0 = (n >> (targetQubit + 1));
    int *ketIndex = new int[n];
    for (i = 0; i < n; i++) {
        ketIndex[i] = i ^ tmp0;
    }
    FastPauli::apply_permutation(ket, ketIndex, n);
    delete[] ketIndex;
}

void FastPauli::ketZ_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit) {
    int i, ketIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    for (i = 0; i < n2; i++) {
        ketIndex = n + 1 - (2 * (i + 1) - tmp0 + FastPauli::fastMod(-i - 1, numQubits - targetQubit - 1)) - 1;
        ket[ketIndex] = -ket[ketIndex];
    }
}

void FastPauli::ketY_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit) {  // Up to some global phase!
    FastPauli::ketZ_r(ket, numQubits, targetQubit);
    FastPauli::ketX_r(ket, numQubits, targetQubit);
}

void FastPauli::ketH_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit) {
    const int n = 1 << numQubits;  // Total size of state vector (2^numQubits)
    // For leftmost indexing, shift by (numQubits - 1 - targetQubit)
    const int targetPower = 1 << (numQubits - 1 - targetQubit);  
    const dtype normFactor = 1.0 / std::sqrt(2.0);  // 1/√2 normalization factor
    
    // Process pairs of elements that differ only in targetQubit
    for (int i = 0; i < n; i++) {
        // Skip if targetQubit bit is set to ensure we process each pair only once
        if (i & targetPower) continue;
        
        // j is i with targetQubit bit flipped
        const int j = i | targetPower;
        
        // Store original values
        const dtype vi = ket[i];
        const dtype vj = ket[j];
        
        // Apply Hadamard transformation:
        // |0⟩ -> (|0⟩ + |1⟩)/√2
        // |1⟩ -> (|0⟩ - |1⟩)/√2
        ket[i] = (vi + vj) * normFactor;
        ket[j] = (vi - vj) * normFactor;
    }
}

void FastPauli::X_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit) {
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!
    int i;
    int j;

    int col;
    int row;
    int tmp;
    int *cIndex = new int[cArrLen];
    const int n = 1 << numQubits;  // Matrix dimension.
    const int tmp0 = (n >> (targetQubit + 1));
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cIndex) private(i, j, col, row, tmp) schedule(static, std::max(1, n / NUMTHREADS))
    // #endif
    for (i = 0; i < n; i++) {
        col = i ^ tmp0;
        for (j = i + 1; j <= n; j++) {
            row = (j - 1) ^ tmp0;
            tmp = FastPauli::getRow(i, j - 1, n);
            cIndex[tmp] = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            // counter++;
            // cArr[cIndex] = -cArr[cIndex];
            // std::cout << "\033[1;33m(col,row): (" << col << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
    FastPauli::apply_permutation(cArr, cIndex, cArrLen);
    delete[] cIndex;
}

void FastPauli::Z_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit) {
    int i;
    int j;
    int col;
    int row;
    int cIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArr) private(i, j, col, row, cIndex) schedule(static, std::max(1, n2 / NUMTHREADS))
    // #endif
    for (i = 0; i < n2; i++) {
        col = 2 * (i + 1) - tmp0 + FastPauli::fastMod(-i - 1, numQubits - targetQubit - 1) - 1;
        // std::cout << "col: " << col << "\n"; // Debug statement.
        for (j = 0; j < n2; j++) {
            row = n + 1 - (2 * (j + 1) - tmp0 + FastPauli::fastMod(-j - 1, numQubits - targetQubit - 1)) - 1;
            cIndex = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            cArr[cIndex] = -cArr[cIndex];
            // std::cout << "\033[1;33m(col,row): (" << col << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
}

void FastPauli::Y_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit) {
    FastPauli::Z_r(cArr, numQubits, targetQubit);
    FastPauli::X_r(cArr, numQubits, targetQubit);
}

// void FastPauli::Y_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit) {
//     const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));
//     int i, col, row, tmp;
//     int *cIndex = new int[cArrLen];
//     const int n = 1 << numQubits;  // Matrix dimension.
//     const int targetBit = 1 << (numQubits - 1 - targetQubit);
    
//     // Build permutation table for Y transformation
//     // Y flips target qubit (like X does) - this is the X part of Y = XZ
//     for (i = 0; i < n; i++) {
//         col = i ^ targetBit;  // Flip target bit
        
//         for (int j = i + 1; j <= n; j++) {  // j is 1-indexed in the loop
//             row = (j - 1) ^ targetBit;  // Flip target bit
            
//             tmp = FastPauli::getRow(i, j - 1, n);
//             cIndex[tmp] = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
//         }
//     }
    
//     // Apply permutation for X part
//     FastPauli::apply_permutation(cArr, cIndex, cArrLen);
    
//     // Apply sign flips from Z part (Y = XZ, global phase i ignored)
//     const int n2 = n >> 1;
//     const int tmp0 = (n >> (targetQubit + 1));
    
//     for (i = 0; i < n2; i++) {
//         col = 2 * (i + 1) - tmp0 + FastPauli::fastMod(-i - 1, numQubits - targetQubit - 1) - 1;
//         for (int j = 0; j < n2; j++) {
//             row = n + 1 - (2 * (j + 1) - tmp0 + FastPauli::fastMod(-j - 1, numQubits - targetQubit - 1)) - 1;
//             int cIndex_val = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
//             cArr[cIndex_val] = -cArr[cIndex_val];
//         }
//     }
    
//     delete[] cIndex;
// }

void FastPauli::H_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit) {
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));  // numQubits cannot be too high, otherwise overflow and cArrLen becomes negative!
    int ii, iii, jj, tmpInd;
    int tmpcounter;
    int col;
    int row;
    int cIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    std::vector<dtype> zz(cArrLen);
    std::vector<dtype> zx(cArrLen);
    std::vector<dtype> xz(cArrLen);
    
    zz = cArr;
    
    // #if defined(_OPENMP)
    //     #pragma omp simd parallel for shared(xz, cArr) private(col, tmpcounter, cIndex) schedule(static)
    // #endif
    for (int i = 0; i < n; i++) {
        col = i ^ tmp0;
        // std::cout << "col: " << col << "\n";
        for (int j = i; j < n; j++) {
            cIndex = (j > col) ? FastPauli::getRow(col, j, n) : FastPauli::getRow(j, col, n);
            tmpcounter = FastPauli::getRow(i, j, n);
            xz[tmpcounter] = cArr[cIndex];
            // counter++;
            // std::cout << "\033[1;33m(col,row): (" << col << "," << j << ")\033[0m\n"; // Debug statement.
        }
    }
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(xz) private(i, ii, j, iii, jj, tmpInd) schedule(static, std::max(1, n2 / NUMTHREADS))
    // #endif
    for (int i = 1; i <= n2; i++) {
        ii = FastPauli::congMod(i - 1, numQubits - targetQubit - 1) + 1;
        for (int j = ii; j <= n; j++) {
            iii = n + 1 - ii;
            jj = n + 1 - j;
            tmpInd = FastPauli::getRow(jj - 1, iii - 1, n);
            xz[tmpInd] = -xz[tmpInd];
            // std::cout << "\033[1;33m(jj,iii): (" << jj-1 << "," << iii-1 << ")\033[0m\n"; // Debug statement.
        }
    }

    // counter = 0;
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(zx, cArr) private(i, j, row, tmpcounter, cIndex) schedule(static, std::max(1, n / NUMTHREADS))
    // #endif
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            row = j ^ tmp0;
            cIndex = (row > i) ? FastPauli::getRow(i, row, n) : FastPauli::getRow(row, i, n);
            tmpcounter = FastPauli::getRow(i, j, n);
            zx[tmpcounter] = cArr[cIndex];
            // counter++;
            // std::cout << "\033[1;33m(col,row): (" << i << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(zx) private(i, ii, j, tmpInd) schedule(static, std::max(1, n2 / NUMTHREADS))
    // #endif
    for (int i = 1; i <= n2; i++) {
        ii = n - FastPauli::congMod(i - 1, numQubits - targetQubit - 1);
        for (int j = ii; j <= n; j++) {
            tmpInd = FastPauli::getRow(ii - 1, j - 1, n);
            zx[tmpInd] = -zx[tmpInd];
            // std::cout << "\033[1;33m(ii,j): (" << ii-1 << "," << j-1 << ")\033[0m\n"; // Debug statement.
        }
    }

    FastPauli::X_r(cArr, numQubits, targetQubit);  // Use original cArr as xx to reduce std::copy usage.
    FastPauli::Z_r(zz, numQubits, targetQubit);
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArr, xz, zx, zz) private(i) schedule(static, std::max(1, cArrLen / NUMTHREADS))
    // #endif
    for (int i = 0; i < cArrLen; i++) {
        cArr[i] = (cArr[i] + xz[i] + zx[i] + zz[i]) / 2.0;
    }

}

std::vector<dtype> FastPauli::Reset_r(const std::vector<dtype>& cArr, const int numQubits, const int targetQubit) {
    std::vector<dtype> outcArr(FastPauli::getCArrLen(numQubits), 0);  // Initialise to zero.
    const int n = 1 << (numQubits - 1);
    const int n2 = 1 << (numQubits - targetQubit - 1);
    int tmpI, tmpJ, tmpIndex, tmpIndex2;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            tmpI = FastPauli::congMod(i, numQubits - targetQubit - 1);
            tmpJ = FastPauli::congMod(j, numQubits - targetQubit - 1);
            tmpIndex = FastPauli::getRow(tmpI, tmpJ, n << 1);
            tmpIndex2 = FastPauli::getRow(tmpI + n2, tmpJ + n2, n << 1);
            outcArr[tmpIndex] = cArr[tmpIndex] + cArr[tmpIndex2];
        }
    }
    return outcArr;
}

void FastPauli::ketCZ_r(std::vector<dtype>& ket, const int numQubits, const int controlQubit, const int targetQubit) {
    const int n = 1 << numQubits;  // Total size of state vector (2^numQubits)
    
    // Calculate bit positions (using leftmost indexing convention)
    const int controlPower = 1 << (numQubits - 1 - controlQubit);
    const int targetPower = 1 << (numQubits - 1 - targetQubit);
    
    // Iterate through all basis states
    for (int i = 0; i < n; i++) {
        // Check if both control and target qubits are 1
        if ((i & controlPower) && (i & targetPower)) {
            // Apply phase flip (-1) when both qubits are 1
            ket[i] = -ket[i];
        }
    }
}

void FastPauli::ketCX_r(std::vector<dtype>& ket, const int numQubits, const int controlQubit, const int targetQubit) {
    const int n = 1 << numQubits;  // Total size of state vector (2^numQubits)
    
    // Calculate bit positions (using leftmost indexing convention)
    const int controlPower = 1 << (numQubits - 1 - controlQubit);
    const int targetPower = 1 << (numQubits - 1 - targetQubit);
    
    // Iterate through all basis states
    for (int i = 0; i < n; i++) {
        // Check if control qubit is 1 and we haven't processed this pair yet
        if ((i & controlPower) && !(i & targetPower)) {
            // j is i with targetQubit flipped
            const int j = i | targetPower;
            
            // Swap amplitudes when control is 1
            std::swap(ket[i], ket[j]);
        }
    }
}

void FastPauli::CZ_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit) {
    assert(controlQubit < targetQubit);
    int i;
    int j;
    int jj;
    int col;
    int row;
    int tmp3;
    int tmp4;
    int tmp5;
    int tmp6;
    int tmp7;
    int tmp8;
    int cIndex;
    const int distance = targetQubit - controlQubit;
    const int totalWidth = 1 << (numQubits - 2);
    const int totalHeight = 3 * totalWidth;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n1 = n >> 1;
    const int tmp0 = n1 >> distance;  // n1/(2^distance).
    const int tmp1 = n >> distance;   // n/(2^distance).

    const int h1 = n1 + tmp0;  // Height of each divided pillar.
    const int h1p = n1 + tmp1;
    const int h2 = totalHeight - h1;  // Height of the checker part (excluding white spaces).

    // std::cout << "h1: " << h1 << ", h1p: " << h1p << ", h2: " << h2 << "\n"; // Debug.

    tmp3 = 1 + controlQubit + targetQubit - numQubits;
    tmp8 = (tmp3 >= 0) ? tmp3 : -tmp3;
    tmp4 = 1 << tmp8;
    tmp5 = (1 << controlQubit) + (1 << targetQubit);
    tmp7 = n >> (targetQubit + 1);
    // std::cout << "tmp3: " << tmp3 << ", tmp8: " << tmp8 << ", tmp4: " << tmp4 << ", tmp5: " << tmp5 << ", tmp7: " << tmp7 << "\n";
    // OPENMP https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArr) private(i, j, col, tmp6, row, cIndex) schedule(static, std::max(1, totalWidth / NUMTHREADS))
    // #endif
    for (i = 0; i < totalWidth; i++) {  // Deal with pillar section of the matrix.
        col = n - FastPauli::congMod(FastPauli::congMod(i, numQubits - targetQubit - 1), numQubits - controlQubit - 1) - 1;
        // std::cout << "col: " << col << "\n"; // Debug statement.
        for (j = 0; j < h1; j++) {
            tmp6 = (tmp3 >= 0) ? (tmp4 * (j)) / tmp5 : (j) / (tmp4 * tmp5);
            row = j + tmp7 * (-1 + (1 << distance)) * tmp6;
            // if (i == totalWidth-1){
            //     std::cout << "row: " << row << "\n"; // Debug statement.
            // }
            cIndex = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            cArr[cIndex] = -cArr[cIndex];
        }
    }
    // #if defined(_OPENMP)
    //     #pragma omp parallel for shared(cArr) private(i, j, col, row, cIndex, jj) schedule(static, std::max(1, totalWidth / NUMTHREADS))
    // #endif
    for (i = 0; i < totalWidth; i++) {  // Deal with checker section of the matrix.
        col = n - FastPauli::congMod(FastPauli::congMod(i, numQubits - targetQubit - 1), numQubits - controlQubit - 1) - 1;
        for (j = h1 + 1; j <= totalHeight; j++) {
            jj = j - h1;
            row = -(tmp7) + 2 * (jj) + FastPauli::fastMod(-(jj), numQubits - targetQubit - 1) + (1 + ((1 << (controlQubit)) * (jj - 1)) / h2) * h1p / (1 << (controlQubit)) - 1;
            cIndex = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
            cArr[cIndex] = -cArr[cIndex];
            // std::cout << "\033[1;33m(col,row): (" << col << "," << row << ")\033[0m\n"; // Debug statement.
        }
    }
}

void FastPauli::CX_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit) {
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));
    int i, col, row, tmp;
    int *cIndex = new int[cArrLen];
    const int n = 1 << numQubits;  // Matrix dimension.
    const int controlBit = 1 << (numQubits - 1 - controlQubit);
    const int targetBit = 1 << (numQubits - 1 - targetQubit);
    
    // Build permutation table for CX transformation
    // CX flips target qubit when control = 1
    for (i = 0; i < n; i++) {
        col = i;
        if (i & controlBit) col ^= targetBit;  // If control=1, flip target bit
        
        for (int j = i + 1; j <= n; j++) {  // j is 1-indexed in the loop
            row = (j - 1);
            if ((j - 1) & controlBit) row ^= targetBit;  // If control=1, flip target bit
            
            tmp = FastPauli::getRow(i, j - 1, n);
            cIndex[tmp] = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
        }
    }
    
    FastPauli::apply_permutation(cArr, cIndex, cArrLen);
    delete[] cIndex;
}

void FastPauli::CY_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit) {
    const int cArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));
    int i, j, col, row, tmp;
    int *cIndex = new int[cArrLen];
    const int n = 1 << numQubits;  // Matrix dimension.
    const int controlBit = 1 << (numQubits - 1 - controlQubit);
    const int targetBit = 1 << (numQubits - 1 - targetQubit);
    
    // Build permutation and tracking for CY transformation
    // CY applies Y to target when control = 1, which swaps |0⟩↔|1⟩ on target
    for (i = 0; i < n; i++) {
        col = i;
        if (i & controlBit) col ^= targetBit;  // If control=1, flip target bit
        
        for (j = i + 1; j <= n; j++) {  // j is 1-indexed in the loop
            row = (j - 1);
            if ((j - 1) & controlBit) row ^= targetBit;  // If control=1, flip target bit
            
            tmp = FastPauli::getRow(i, j - 1, n);
            cIndex[tmp] = (row > col) ? FastPauli::getRow(col, row, n) : FastPauli::getRow(row, col, n);
        }
    }
    
    // Apply permutation for the basis state transformation
    FastPauli::apply_permutation(cArr, cIndex, cArrLen);
    
    // Apply sign flips from Y = iXZ (ignoring global phase i)
    // When target flips for control=1 states in off-diagonal elements, apply sign flip
    for (i = 0; i < n; i++) {
        if (!(i & controlBit)) continue;  // Skip control=0 states
        
        for (j = i + 1; j < n; j++) {
            // Check if target bits differ in the original pair (i, j)
            int i_target = (i & targetBit);
            int j_target = (j & targetBit);
            
            // If target differs and control=1 for both, apply sign flip to off-diagonal
            if (i_target != j_target && (j & controlBit)) {
                int cIndex_val = FastPauli::getRow(i, j, n);
                cArr[cIndex_val] = -cArr[cIndex_val];
            }
        }
    }
    
    delete[] cIndex;
}

// Caculates trace of given DM (in the form of cArr)
dtype FastPauli::trace_r(const std::vector<dtype>& cArr, const int numQubits) {
    dtype out = 0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        out += cArr[FastPauli::getRow(i, i, n)];
    }
    return out;
}

dtype FastPauli::trace_ket(const std::vector<dtype>& ket){
    dtype out = 0.0;
    for (auto const &el : ket)
    {
        out += el * el;
    }
    return out;
}

dtype FastPauli::trace2_r(const std::vector<dtype>& cArr, const int numQubits) {
    dtype out = 0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            out += 2 * cArr[FastPauli::getRow(i, j, n)] * cArr[FastPauli::getRow(i, j, n)];
        }
    }
    for (int i = 0; i < n; i++) {
        out += cArr[FastPauli::getRow(i, i, n)] * cArr[FastPauli::getRow(i, i, n)];
    }
    return out;
}

dtype FastPauli::fidel_r(const std::vector<dtype>& ket, const std::vector<dtype>& cArr, const int numQubits) {
    dtype out = 0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            out += 2 * ket[j] * ket[i] * cArr[FastPauli::getRow(i, j, n)];
        }
    }
    for (int i = 0; i < n; i++) {
        out += ket[i] * ket[i] * cArr[FastPauli::getRow(i, i, n)];
    }
    return out;
}

dtype FastPauli::fidel_ket(const std::vector<dtype>& in_ket, const std::vector<dtype>& corr_ket, const int numQubits) {
    dtype bra_ket = 0.0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        bra_ket += in_ket[i] * corr_ket[i];
    }
    return bra_ket * bra_ket;
}

void FastPauli::printUTri(const std::vector<dtype>& cArr, const int numQubits) {
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j >= i) {
                std::cout << "\033[1m";
                std::cout << std::setw(5) << cArr[FastPauli::getRow(i, j, n)];
                std::cout << "\033[0m";
                std::cout << " ";
            } else {
                std::cout << std::setw(5) << 0;
                std::cout << " ";
            }
        }
        std::cout << "\n";
    }
}