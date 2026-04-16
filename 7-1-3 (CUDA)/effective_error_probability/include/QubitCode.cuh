// Changes related to number of encoded qubits and STABILIZERS must be made
// '16' NEEDS TO BE CHANGED

#pragma once

#include <algorithm> // Provide std::transform. Reference: https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint> // Provides definition for uint64_t.
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <tuple>
#include <format>

#include "prereq.cuh"
#include "fpauli.cuh"
#include "npauli.cuh"

namespace QubitCode
{

    /**
     * @brief Applies correction to the input density matrix according to the syndrome `i`. Not fault-tolerant.
     *
     * @param i i ∈ {0,1,2,...,DIM_ANCILLA-1} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     */
    void applyCorrection(const uint16_t i, std::vector<double> &cArr); // For non-fault tolerant implementation of the [[7,1,3]] code, i.e., no flag qubits.

    /**
     * @brief Applies correction to the input density matrix according to the syndrome `i` and which `flag ∈ {0,1,2,3,4,5}` was triggered. Fault-tolerant.
     * See Bernard's MSc thesis.
     *
     * @param i i ∈ {0,1,2,...,DIM_ANCILLA-1} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param flag flag ∈ {0,1,2,3,4,5} specifies which flag was triggered.
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     */
    void applyCorrection_flag(const uint16_t i, const uint16_t flag, std::vector<double> &cArr);

    /**
     * @brief Applies 1-erasure error correction to the input density matrix according to the syndrome `i`. Not fault-tolerant by nature.
     *
     * @param i i ∈ {0,1,2,...,63} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     * @param k i ∈ {0,1,2,3,4,5,6} is qubit lost index.
     * DEFINED in ..._loss.cpp
     */
    void applyCorrection_1erasure(const uint16_t i, std::vector<double> &cArr, const int k);

    enum class correctionMode
    {
        noLossW1,
        noLossW2,
        twoLoss,
        oneLoss
    };

    struct ECInfo
    {
        std::vector<double> grandErrors;
        std::vector<double> grandTraces;
        void resize(size_t numElems)
        {
            grandErrors.resize(numElems);
            grandTraces.resize(numElems);
        }
        void save(std::string const &filename)
        {
            std::ofstream outfile_noappend;
            std::ofstream outfile_append;
            outfile_noappend.open(filename, std::ios::binary);
            outfile_noappend.write((char *)grandErrors.data(), grandErrors.size() * sizeof(double));
            outfile_noappend.close();

            outfile_append.open(filename, std::ios::binary | std::ios::app);
            outfile_append.write((char *)grandTraces.data(), grandTraces.size() * sizeof(double));
            outfile_append.close();
        }
        void print()
        {
            std::cout << "\x1b[0;95mecInfo.grandTraces = {";
            for (size_t i = 0; i < grandTraces.size(); i++)
            {
                if (i == grandTraces.size() - 1)
                {
                    std::cout << grandTraces[i];
                }
                else
                {
                    std::cout << grandTraces[i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;

            std::cout << "====================================" << std::endl;
            std::cout << "\x1b[0;94mecInfo.grandErrors = {";
            for (size_t i = 0; i < grandErrors.size(); i++)
            {
                if (i == grandErrors.size() - 1)
                {
                    std::cout << grandErrors[i];
                }
                else
                {
                    std::cout << grandErrors[i] << ", ";
                }
            }
            std::cout << "};\x1b[0m" << std::endl;
        }
    };

    Eigen::Matrix2f getPauli(const char p);

    void Form_Satbilizers(const std::string &filename);

    /**
     * @brief The first set of stabiliser operations corresponding to `XZZXI`. Fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)`
     * where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    void S1_flag(std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `IXZZX`. Fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)`
     * where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    void S2_flag(std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `XIXZZ`. Fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)`
     * where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    void S3_flag(std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `ZXIXZ`. Fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)`
     * where `n=7` is the number of qubits (including the ancilla and flag qubit).
     * @param err The error rate of the two-qubit gates in the stabiliser operation between the ancilla and the data qubits, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     * @param ancillaFlagErr The error rate of the two-qubit gates between the ancilla and the flag qubit. Usually same as `err`.
     */
    void S4_flag(std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    // Can add desc. if you want. And above descs are only for 5-1-3
    void S5_flag(std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    void S6_flag(std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    void Sn_flag(int n, std::vector<double> &cArr, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief The first set of stabiliser operations corresponding to `IIIXXXX`. Not fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    void S1Sequential(std::vector<double> &cArr, const double err, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `IXXIIXX`. Not fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    void S2Sequential(std::vector<double> &cArr, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `XIXIXIX`. Not fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    void S3Sequential(std::vector<double> &cArr, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `IIIZZZZ`. Not fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    void S4Sequential(std::vector<double> &cArr, const double errTeleport);

    /**
     * @brief The second set of stabiliser operations corresponding to `IZZIIZZ`. Not fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    void S5Sequential(std::vector<double> &cArr, const double errTeleport);


    /**
     * @brief The second set of stabiliser operations corresponding to `ZIZIZIZ`. Not fault-tolerant.
     *
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=8` is the number of qubits.
     * @param err The error rate of the two-qubit gates in the stabiliser operation, i.e., the CNOT and CPHASE gates.
     * @param errTeleport The error rate of the teleported two-qubit gates in 3 of the 4 two-qubit gates. This is typically higher than `err`.
     */
    void S6Sequential(std::vector<double> &cArr, const double errTeleport);

    void Sn_Sequential(int n, std::vector<double> &cArr, const double err, const double errTeleport);

    // The 3 below ARE FINE AS IS
    void nested_loops(int ancilla, int &corr_ind, std::vector<std::vector<double>> &cArrN, std::vector<double> &cArr_Nplus1, std::vector<double> &result, const std::vector<std::vector<double>> &ketH, const std::vector<double> &cArrH1P, const double err, const double errTeleport, const QubitCode::correctionMode correctionMode, const int qblost_or_flag);

    std::vector<double> OneRoundSequential(const std::vector<double> &cArrCodeword, const double err, const double errTeleport, const QubitCode::correctionMode correctionMode, const int qblost_or_flag);

    std::vector<double> OneRoundSequential_flag(const std::vector<double> &logicalcArr, const std::vector<double> &cArrCodeword, const double err, const double errTeleport, const double ancillaFlagErr);

    /**
     * @brief Returns an array that corresponds to the cArr vector |+⟩^(⊗numQubits).
     *
     * @tparam numQubits
     * @tparam Size
     * @return std::vector<double, Size>
     */
    template <int numQubits>
    std::vector<double> getKetH()
    { // Gets a cArr vector.
        const double elemVal = 1.0 / std::pow(std::sqrt(2), numQubits);
        std::vector<double> Ket(FastPauli::KetLen<numQubits>, elemVal);
        return Ket;
    }

    /**
     * @brief Get the std::vector<double, 10> result which represents
     * DensityMatrix(ancilla,flag). For example, if `isAncillaTriggered = true` and `isFlagTriggered = true`, then the result is DensityMatrix(+,0).
     *
     * @tparam isAncillaTriggered true if ancilla triggered, else false.
     * @tparam isFlagTriggered true if flag triggered, else false.
     * @return constexpr std::vector<double, FastPauli::CArrLen<2>>
     */
    template <bool isAncillaTriggered, bool isFlagTriggered> // NO CHANGE NEEDED. cArr versions of given DMs
    std::vector<double> getAncillaFlagDensityMatrix()
    {
        if (!isAncillaTriggered && !isFlagTriggered)
        { // DensityMatrix(+,0)
            return {0.5, 0, 0.5, 0, 0, 0, 0, 0.5, 0, 0};
        }
        else if (isAncillaTriggered && !isFlagTriggered)
        { // DensityMatrix(-,0)
            return {0.5, 0, -0.5, 0, 0, 0, 0, 0.5, 0, 0};
        }
        else if (!isAncillaTriggered && isFlagTriggered)
        { // DensityMatrix(+,1)
            return {0, 0, 0, 0, 0.5, 0, 0.5, 0, 0, 0.5};
        }
        else
        { // DensityMatrix(-,1)
            return {0, 0, 0, 0, 0.5, 0, -0.5, 0, 0, 0.5};
        }
    }

    template <bool isAncillaTriggered, bool isFlagTriggered>
    std::vector<double> getAncillaFlagKet() // NO CHANGE NEEDED. Ket version of previous function.
    {
        if (!isAncillaTriggered && !isFlagTriggered)
        { // Ket(+,0)
            return {std::sqrt(2.0) / 2, 0, std::sqrt(2.0) / 2, 0};
        }
        else if (isAncillaTriggered && !isFlagTriggered)
        { // Ket(-,0)
            return {std::sqrt(2.0) / 2, 0, -std::sqrt(2.0) / 2, 0};
        }
        else if (!isAncillaTriggered && isFlagTriggered)
        { // Ket(+,1)
            return {0, std::sqrt(2.0) / 2, 0, std::sqrt(2.0) / 2};
        }
        else
        { // Ket(-,1)
            return {0, std::sqrt(2.0) / 2, 0, -std::sqrt(2.0) / 2};
        }
    }

    /**
     * @brief Get the Logical Ket vector for the [[7,1,3]] code given by |ψ_L⟩ = cos(θ/2)|0_L⟩ + sin(θ/2)|1_L⟩.
     *
     * @param theta A real number θ such that |ψ_L⟩ = cos(θ/2)|0_L⟩ + sin(θ/2)|1_L⟩.
     * @return std::vector<double, 128>
     */
    std::vector<datatype> getLogicalKet(const double theta);

    std::vector<datatype> getLogicalDensityMatrix(const double theta);

    /**
     * @brief This function basically does (⟨ancilla|⊗𝕀...⊗𝕀).ρ.(|ancilla⟩⊗𝕀...⊗𝕀)
     * where ρ is the input real density matrix and 𝕀 is a 2×2 identity matrix.
     * Effectively, what this does is gives you the cArr for data qubits (extracted from ρ) ONLY when the input ancilla MATCHES those of ρ
     * Otherwise the result will become 0
     *
     * @param in_cArr The output density matrix with reduced dimension of length `pow(2,numQubits-1)+pow(2,2numQubits-1)`.
     * @param numAncillaQubits The number of ancilla qubits.
     * @param numQubits The number of qubits in the resulting reduced density matrix.
     * @param ancillaKet The cArr vector of the ancilla qubits, i.e., |ancilla⟩.
     * @param out_cArr The input density matrix with the full dimension of length
     * `pow(2,(numQubits+numAncillaQubits)-1)+pow(2,2(numQubits+numAncillaQubits)-1)`.
     */
    void applyAncillaKet(std::vector<double> &out_cArr, const int numAncillaQubits, const int numQubits, const std::vector<double> &ancillaKet, const std::vector<double> &in_cArr);

    template <int numAncillaQubits, int numQubits>
    std::vector<double> applyAncillaKet(const std::vector<double> &ancillaKet, const std::vector<double> &in_cArr)
    {
        std::vector<double> out_cArr(FastPauli::CArrLen<numQubits>);
        QubitCode::applyAncillaKet(out_cArr, numAncillaQubits, numQubits, ancillaKet, in_cArr);
        return out_cArr;
    }

    std::vector<double> getcArrW0(const std::vector<double> &cArr);
    std::vector<double> getcArrW1(const std::vector<double> &cArr);
    std::vector<double> getcArrW2(const std::vector<double> &cArr);
};