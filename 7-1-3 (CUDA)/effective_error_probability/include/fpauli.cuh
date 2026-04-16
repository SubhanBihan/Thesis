// Nothing to change here

#pragma once

#include <algorithm>  // Provides std::copy
#include <array>
#include <cassert>
#include <cmath>  // Provides std::max
#include <cstddef>
#include <cstdint>  // Provides definition for uint64_t.
#include <iomanip>
#include <iostream>
#include <vector>
// #include <iterator> // std::begin() and std::end()

#include "AuxMath.cuh"

namespace FastPauli {

    void apply_permutation(std::vector<double>& v, int *indices, const int vSize);
    int fastMod(const int i, const int t);

    template <unsigned int numQubits>
    constexpr uint64_t CArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));

	template <unsigned int numQubits>
    constexpr uint64_t KetLen = (UINT64_C(1) << numQubits);
        

    /**
     * @brief Numbers congruent to mod 2^t. In Matlab, it would be '2*i - 2^t + mod(-i,2^t)'.
     * Note that this uses index starting from 1 as opposed to the C/C++ function here which uses index starting from 0.
     * congMod(i,t) gives the numbers that are congruent to {0,1,...,2^t-1} mod 2^(t+1).
     * @param i Integer.
     * @param t Exponent.
     * @return int
     */
    int congMod(const int i, const int t);
    int getRow(int i, int j, const int n);
    const uint64_t getCArrLen(const int numQubits);

    /**
     * @brief This function basically does `DensityMatrixNew = DensityMatrixOld1 ⊗ DensityMatrixOld2`
     */
    void dmKron(std::vector<double>& result, const std::vector<double>& CArr1, const std::vector<double>& CArr2, int numQubits1, int numQubits2);

    // Kronecker product function that takes and returns vectors
    template <int numQubits1, int numQubits2>
    std::vector<double> dmKron(const std::vector<double>& CArr1, const std::vector<double>& CArr2) {
        std::vector<double> resultvec(Operator_tnumQubits1 + numQubits2>);
        dmKron(resultvec, CArr1, CArr2, numQubits1, numQubits2);

        return resultvec;
    }

    void ketKron(std::vector<double>& result, const std::vector<double>& ket1, const std::vector<double>& ket2, const int numQubits1, const int numQubits2);

    template <int numQubits1, int numQubits2>
    std::vector<double> ketKron(const std::vector<double>& ket1, const std::vector<double>& ket2) {
        std::vector<double> resultvec(KetLen<numQubits1 + numQubits2>);
        ketKron(resultvec, ket1, ket2, numQubits1, numQubits2);

        return resultvec;
    }
    

    void ketToCArr(std::vector<double>& cArr, const std::vector<double>& ket, const int numQubits);
    void ketX_r(std::vector<double>& ket, const int numQubits, const int targetQubit);
    void ketZ_r(std::vector<double>& ket, const int numQubits, const int targetQubit);
    void ketY_r(std::vector<double>& ket, const int numQubits, const int targetQubit);  // Up to some global phase!
    void ketH_r(std::vector<double>& ket, const int numQubits, const int targetQubit);
    void X_r(std::vector<double>& cArr, const int numQubits, const int targetQubit);
    void Z_r(std::vector<double>& cArr, const int numQubits, const int targetQubit);
    void Y_r(std::vector<double>& cArr, const int numQubits, const int targetQubit);
    void H_r(std::vector<double>& cArr, const int numQubits, const int targetQubit);
    std::vector<double> Reset_r(const std::vector<double>& cArr, const int numQubits, const int targetQubit);
    
    template <int numQubits>
    std::vector<double> Reset_r(const std::vector<double>& cArr, const int targetQubit) {
        std::vector<double> outcArr(CArrLen<numQubits>, 0);  // Initialise to zero.
        const int n = 1 << (numQubits - 1);
        const int n2 = 1 << (numQubits - targetQubit - 1);
        int tmpI, tmpJ, tmpIndex, tmpIndex2;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                tmpI = congMod(i, numQubits - targetQubit - 1);
                tmpJ = congMod(j, numQubits - targetQubit - 1);
                tmpIndex = getRow(tmpI, tmpJ, n << 1);
                tmpIndex2 = getRow(tmpI + n2, tmpJ + n2, n << 1);
                outcArr[tmpIndex] = cArr[tmpIndex] + cArr[tmpIndex2];
            }
        }
        return outcArr;
    }

    void CZ_r(std::vector<double>& cArr, const int numQubits, const int controlQubit, const int targetQubit);
    void CX_r(std::vector<double>& cArr, const int numQubits, const int controlQubit, const int targetQubit);
    void ketCZ_r(std::vector<double>& ket, const int numQubits, const int controlQubit, const int targetQubit);
    void ketCX_r(std::vector<double>& ket, const int numQubits, const int controlQubit, const int targetQubit);
    double fidel_r(const std::vector<double>& ket, const std::vector<double>& cArr, const int numQubits);
    double fidel_ket(const std::vector<double>& in_ket, const std::vector<double>& corr_ket, const int numQubits);
    double trace_r(const std::vector<double>& cArr, const int numQubits);
    double trace_ket(const std::vector<double>& ket);
    double trace2_r(const std::vector<double>& cArr, const int numQubits);
    void printUTri(const std::vector<double>& cArr, const int numQubits);

};