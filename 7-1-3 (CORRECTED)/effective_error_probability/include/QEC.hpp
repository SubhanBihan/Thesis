// Changes related to number of encoded qubits and STABILIZERS must be made
// '16' NEEDS TO BE CHANGED

#ifndef _FIVEQUBITCODEHPP_
#define _FIVEQUBITCODEHPP_

#include <algorithm> // Provide std::transform. Reference: https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
#include <array>
#include <cstddef>
#include <cstdint> // Provides definition for uint64_t.
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>
#include <format>

#include "params.hpp"
#include "fpauli.hpp"
#include "npauli.hpp"

extern std::vector<dtype> ketH1P;

// Get ready |+⟩⟨+|.
extern std::vector<dtype> cArrH1P;

// Get ready |-⟩.
extern std::vector<dtype> ketH1N;

extern std::vector<std::vector<dtype>> ketH;

class QEC
{
private:
    /**
     * @brief Applies correction to the input density matrix according to the syndrome `i`. Not fault-tolerant.
     *
     * @param i i ∈ {0,1,2,...,63} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     */
    static void applyCorrection(const uint16_t i, std::vector<dtype> &cArr); // For non-fault tolerant implementation of the [[7,1,3]] code, i.e., no flag qubits.

    /**
     * @brief Applies correction to the input density matrix according to the syndrome `i` and which `flag ∈ {0,1,2,3,4,5}` was triggered. Fault-tolerant.
     * See Bernard's MSc thesis.
     *
     * @param i i ∈ {0,1,2,...,63} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param flag flag ∈ {0,1,2,3,4,5} specifies which flag was triggered.
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     */
    static void applyCorrection_flag(const uint16_t i, const uint16_t flag, std::vector<dtype> &cArr);

    /**
     * @brief Applies 1-erasure error correction to the input density matrix according to the syndrome `i`. Not fault-tolerant by nature.
     *
     * @param i i ∈ {0,1,2,...,63} is all the possibilities of the ancilla qubits on/off combinations, i.e., the syndromes.
     * @param cArr The real density matrix in the form of a C-style array of length `pow(2,n-1)+pow(2,2n-1)` where `n=7` is the number of qubits.
     * @param k i ∈ {0,1,2,3,4,5,6} is qubit lost index.
     * DEFINED in ..._loss.cpp
     */
    static void applyCorrection_1erasure(const uint16_t i, std::vector<dtype> &cArr, const int k);

public:
    enum correctionMode // NOT SURE whether any changes/additions are necessary
    {
        noLossW1,
        noLossW2,
        oneLoss
    };

    struct ECInfo
    {
        // std::array<std::vector<dtype>, NUMQUBITS> errors; // 1 - fidelity.
        // std::array<std::vector<dtype>, NUMQUBITS> traces; // These act like probabilities.
        std::vector<dtype> grandErrors;
        std::vector<dtype> grandTraces;
        void resize(size_t numElems)
        {
            // for (int i = 0; i < NUMQUBITS; i++)
            // {
            //     errors[i].resize(numElems);
            //     traces[i].resize(numElems);
            // }
            grandErrors.resize(numElems);
            grandTraces.resize(numElems);
        }
        void save(std::string const &filename)
        {
            std::ofstream outfile_noappend;
            std::ofstream outfile_append;
            outfile_noappend.open(filename, std::ios::binary);
            outfile_noappend.write((char *)grandErrors.data(), grandErrors.size() * sizeof(dtype));
            outfile_noappend.close();

            outfile_append.open(filename, std::ios::binary | std::ios::app);
            outfile_append.write((char *)grandTraces.data(), grandTraces.size() * sizeof(dtype));
            outfile_append.close();
            
            // outfile_noappend.open(filename, std::ios::binary);
            // outfile_noappend.write((char *)errors[0].data(), errors[0].size() * sizeof(dtype));
            // outfile_noappend.close();

            // outfile_append.open(filename, std::ios::binary | std::ios::app);
            // for (int i = 1; i < NUMQUBITS; i++)
            // {
            //     outfile_append.write((char *)errors[i].data(), errors[i].size() * sizeof(dtype));
            // }
            // for (int i = 0; i < NUMQUBITS; i++)
            // {
            //     outfile_append.write((char *)traces[i].data(), traces[i].size() * sizeof(dtype));
            // }
            // outfile_append.write((char *)grandErrors.data(), grandErrors.size() * sizeof(dtype));
            // outfile_append.write((char *)grandTraces.data(), grandTraces.size() * sizeof(dtype));
            // outfile_append.close();
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

            // for (size_t j = 0; j < NUMQUBITS; j++)
            // {
            //     std::cout << std::format("\x1b[0;92mecInfo.traces[{}] = {{", j);
            //     for (size_t i = 0; i < traces[j].size(); i++)
            //     {
            //         if (i == traces[j].size() - 1)
            //         {
            //             std::cout << traces[j][i];
            //         }
            //         else
            //         {
            //             std::cout << traces[j][i] << ", ";
            //         }
            //     }
            //     std::cout << "};\x1b[0m" << std::endl;
            // }

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

            // for (size_t j = 0; j < NUMQUBITS; j++)
            // {
            //     std::cout << std::format("\x1b[0;92mecInfo.errors[{}] = {{", j);
            //     for (size_t i = 0; i < errors[j].size(); i++)
            //     {
            //         if (i == errors[j].size() - 1)
            //         {
            //             std::cout << errors[j][i];
            //         }
            //         else
            //         {
            //             std::cout << errors[j][i] << ", ";
            //         }
            //     }
            //     std::cout << "};\x1b[0m" << std::endl;
            // }
        }
    };

    // Flagged
    static void S1_flag(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    static void S2_flag(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    static void S3_flag(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    static void S4_flag(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    static void S5_flag(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    static void S6_flag(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    static void Sn_flag(int n, std::vector<dtype> &cArr, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);

    // Non-FT
    static void S1Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);

    static void S2Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);

    static void S3Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);

    static void S4Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);

    static void S5Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);

    static void S6Sequential(std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);

    static void Sn_Sequential(int n, std::vector<dtype> &cArr, const dtype err, const dtype errTeleport);
   

    // Critical APIs
    static void ORS_NestedLoops(int ancilla, int &corr_ind, std::vector<std::vector<dtype>> &cArrN, std::vector<dtype> &result, const dtype err, const dtype errTeleport, const QEC::correctionMode correctionMode, const int qblost_or_flag);

    static std::vector<dtype> OneRoundSequential(const std::vector<dtype> &cArrCodeword, const dtype err, const dtype errTeleport, const QEC::correctionMode correctionMode, const int qblost_or_flag);

    static std::vector<dtype> OneRoundSequential_flag(const std::vector<dtype> &logicalcArr, const std::vector<dtype> &cArrCodeword, const dtype err, const dtype errTeleport, const dtype ancillaFlagErr);


    /**
     * @brief Returns an array that corresponds to the cArr vector |+⟩^(⊗numQubits).
     *
     * @tparam numQubits
     * @tparam Size
     * @return std::vector<dtype, Size>
     */
    template <int numQubits>
    static std::vector<dtype> getKetH()
    { // Gets a cArr vector.
        const dtype elemVal = 1.0 / std::pow(std::sqrt(2), numQubits);
        std::vector<dtype> Ket(FastPauli::KetLen<numQubits>, elemVal);
        return Ket;
    }


    /**
     * @brief Get the std::vector<dtype, 10> result which represents
     * DensityMatrix(ancilla,flag). For example, if `isAncillaTriggered = true` and `isFlagTriggered = true`, then the result is DensityMatrix(+,0).
     *
     * @tparam isAncillaTriggered true if ancilla triggered, else false.
     * @tparam isFlagTriggered true if flag triggered, else false.
     * @return constexpr std::vector<dtype, FastPauli::CArrLen<2>>
     */
    template <bool isAncillaTriggered, bool isFlagTriggered> // NO CHANGE NEEDED. cArr versions of given DMs
    static std::vector<dtype> getAncillaFlagDensityMatrix()
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
    static std::vector<dtype> getAncillaFlagKet() // NO CHANGE NEEDED. Ket version of previous function.
    {
        if (!isAncillaTriggered && !isFlagTriggered)
        { // Ket(+,0)
            return {num, 0, num, 0};
        }
        else if (isAncillaTriggered && !isFlagTriggered)
        { // Ket(-,0)
            return {num, 0, -num, 0};
        }
        else if (!isAncillaTriggered && isFlagTriggered)
        { // Ket(+,1)
            return {0, num, 0, num};
        }
        else
        { // Ket(-,1)
            return {0, num, 0, -num};
        }
    }

    /**
     * @brief Get the Logical Ket vector for the [[7,1,3]] code given by |ψ_L⟩ = cos(θ/2)|0_L⟩ + sin(θ/2)|1_L⟩.
     *
     * @param theta A real number θ such that |ψ_L⟩ = cos(θ/2)|0_L⟩ + sin(θ/2)|1_L⟩.
     * @return std::vector<dtype, 128>
     */
    static std::vector<dtype> getLogicalKet(const dtype theta);

    static std::vector<dtype> getLogicalDensityMatrix(const dtype theta);

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
    static void applyAncillaKet(std::vector<dtype> &out_cArr, const int numAncillaQubits, const int numQubits, const std::vector<dtype> &ancillaKet, const std::vector<dtype> &in_cArr);

    template <int numAncillaQubits, int numQubits>
    static std::vector<dtype> applyAncillaKet(const std::vector<dtype> &ancillaKet, const std::vector<dtype> &in_cArr)
    {
        std::vector<dtype> out_cArr(FastPauli::CArrLen<numQubits>);
        QEC::applyAncillaKet(out_cArr, numAncillaQubits, numQubits, ancillaKet, in_cArr);
        return out_cArr;
    }

    static std::vector<dtype> getcArrW0(const std::vector<dtype> &cArr);
    static std::vector<dtype> getcArrW1(const std::vector<dtype> &cArr);
    static std::vector<dtype> getcArrW2(const std::vector<dtype> &cArr);
};

#endif