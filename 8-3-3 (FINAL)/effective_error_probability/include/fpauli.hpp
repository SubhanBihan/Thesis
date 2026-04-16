// Nothing to change here
#ifndef _FPAULIHPP_
#define _FPAULIHPP_

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

#include "AuxMath.hpp"

class FastPauli {
   private:
    static void apply_permutation(std::vector<dtype>& v, int *indices, const int vSize);
    static int fastMod(const int i, const int t);

   public:
    template <int numQubits>
    static constexpr uint64_t CArrLen = (UINT64_C(1) << (numQubits - 1)) + (UINT64_C(1) << (2 * numQubits - 1));

	template <int numQubits>
    static constexpr uint64_t KetLen = (UINT64_C(1) << numQubits);
        

    /**
     * @brief Numbers congruent to mod 2^t. In Matlab, it would be '2*i - 2^t + mod(-i,2^t)'.
     * Note that this uses index starting from 1 as opposed to the C/C++ function here which uses index starting from 0.
     * congMod(i,t) gives the numbers that are congruent to {0,1,...,2^t-1} mod 2^(t+1).
     * @param i Integer.
     * @param t Exponent.
     * @return int
     */
    static int congMod(const int i, const int t);
    static int getRow(int i, int j, const int n);
    static const uint64_t getCArrLen(const int numQubits);

    /**
     * @brief This function basically does `DensityMatrixNew = DensityMatrixOld1 ⊗ DensityMatrixOld2`
     */
    static void dmKron(std::vector<dtype>& result, const std::vector<dtype>& CArr1, const std::vector<dtype>& CArr2, int numQubits1, int numQubits2);

    // Kronecker product function that takes and returns vectors
    template <int numQubits1, int numQubits2>
    static std::vector<dtype> dmKron(const std::vector<dtype>& CArr1, const std::vector<dtype>& CArr2) {
        std::vector<dtype> resultvec(CArrLen<numQubits1 + numQubits2>);
        dmKron(resultvec, CArr1, CArr2, numQubits1, numQubits2);

        return resultvec;
    }

    static void ketKron(std::vector<dtype>& result, const std::vector<dtype>& ket1, const std::vector<dtype>& ket2, const int numQubits1, const int numQubits2);

    template <int numQubits1, int numQubits2>
    static std::vector<dtype> ketKron(const std::vector<dtype>& ket1, const std::vector<dtype>& ket2) {
        std::vector<dtype> resultvec(KetLen<numQubits1 + numQubits2>);
        ketKron(resultvec, ket1, ket2, numQubits1, numQubits2);

        return resultvec;
    }
    

    static void ketToCArr(std::vector<dtype>& cArr, const std::vector<dtype>& ket, const int numQubits);
    static void ketX_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit);
    static void ketZ_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit);
    static void ketY_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit);  // Up to some global phase!
    static void ketH_r(std::vector<dtype>& ket, const int numQubits, const int targetQubit);
    static void X_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit);
    static void Z_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit);
    static void Y_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit);
    static void H_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit);
    static std::vector<dtype> Reset_r(const std::vector<dtype>& cArr, const int numQubits, const int targetQubit);
    
    template <int numQubits>
    static std::vector<dtype> Reset_r(const std::vector<dtype>& cArr, const int targetQubit) {
        std::vector<dtype> outcArr(CArrLen<numQubits>, 0);  // Initialise to zero.
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

    static void CZ_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit);
    static void CX_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit);
    static void CY_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit);
    static void ketCZ_r(std::vector<dtype>& ket, const int numQubits, const int controlQubit, const int targetQubit);
    static void ketCX_r(std::vector<dtype>& ket, const int numQubits, const int controlQubit, const int targetQubit);
    static dtype fidel_r(const std::vector<dtype>& ket, const std::vector<dtype>& cArr, const int numQubits);
    static dtype fidel_ket(const std::vector<dtype>& in_ket, const std::vector<dtype>& corr_ket, const int numQubits);
    static dtype trace_r(const std::vector<dtype>& cArr, const int numQubits);
    static dtype trace_ket(const std::vector<dtype>& ket);
    static dtype trace2_r(const std::vector<dtype>& cArr, const int numQubits);
    static void printUTri(const std::vector<dtype>& cArr, const int numQubits);

};

#endif