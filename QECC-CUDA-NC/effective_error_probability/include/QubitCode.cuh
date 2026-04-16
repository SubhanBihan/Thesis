#pragma once

#include <array>
#include <cstddef>
#include <cstdint> // Provides definition for uint64_t.
#include <functional>
#include <numbers>
#include <tuple>
#include <format>
#include <thread>

#include "prereq.cuh"
#include "kernels.cuh"

namespace QubitCode {
    enum class correctionMode
    {
        noLossW1,
        noLossW2,
        oneLoss,
        twoLoss
    };

    /**
     * @brief Get the Logical Ket/CArr/Full-DM vector for the QECC given by | _L  = cos( /2)|0_L  + sin( /2)|1_L .
     *
     * @param theta A real number   such that | _L  = cos( /2)|0_L  + sin( /2)|1_L .
     */
    std::vector<datatype> getLogicalKet(const double theta);
    std::vector<datatype> getLogicalDensityMatrix(const double theta);

    void Sn_seq(datatype* dummy, datatype* input, unsigned n, const unsigned numQubits, const datatype err, const datatype errTeleport,
                std::vector<cudaStream_t>& streams, datatype** d_temps);

    void Sn_flag(datatype* dummy, datatype* input, unsigned n, const unsigned numQubits, const datatype err, const datatype errTeleport, const datatype ancillaFlagErr,
                 std::vector<cudaStream_t>& streams, datatype** d_temps);

    void applyAncilla(datatype* output, datatype* ancilla_mat, datatype* input, unsigned m, unsigned n, cublasHandle_t handle);

    void applyCorrection(datatype* dummy, datatype* input, unsigned i, dim3 gridSize, cudaStream_t stream);

    void applyCorrection_1erasure(datatype* dummy, datatype* input, unsigned i, unsigned qblost, dim3 gridSize, cudaStream_t stream);
    
    void applyCorrection_flag(datatype* dummy, datatype* input, unsigned i, unsigned flag, dim3 gridSize, cudaStream_t stream);

    void ORS_NestedLoops(datatype* dummy, unsigned num_anc, unsigned &corr_ind, std::vector<datatype*>& ORS_DMs, datatype* DM_Nplus1, datatype** d_temps, std::vector<cudaStream_t>& streams, cublasHandle_t handle, 
                          dim3 gridSize_N, dim3 gridSize_Nplus1, const datatype err, const datatype errTeleport, const QubitCode::correctionMode correctionMode, const unsigned qblost_or_flag);
                         
    void OneRoundSequential(datatype* dummy, datatype* input, std::vector<datatype*>& ORS_DMs, datatype* DM_Nplus1, datatype** d_temps, std::vector<cudaStream_t>& streams, cublasHandle_t handle, 
                            const datatype err, const datatype errTeleport, const QubitCode::correctionMode correctionMode, const unsigned qblost_or_flag);

    void OneRoundSequential_flag(datatype* input, std::vector<std::vector<datatype*>>& ORS_DMs_sets, std::vector<datatype*>& DMs_Nplus2, std::vector<datatype*>& DMs_Nplus1, std::array<datatype**, 2>& d_temps, std::array<std::vector<cudaStream_t>, 2>& streams_sets, std::array<std::vector<cublasHandle_t>, 2>& handles_sets, 
                            const datatype err, const datatype errTeleport, const datatype ancillaFlagErr);

    struct ECInfo
    {
        std::vector<datatype> grandErrors;
        std::vector<datatype> grandTraces;
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
            outfile_noappend.write((char *)grandErrors.data(), grandErrors.size() * sizeof(datatype));
            outfile_noappend.close();

            outfile_append.open(filename, std::ios::binary | std::ios::app);
            outfile_append.write((char *)grandTraces.data(), grandTraces.size() * sizeof(datatype));
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
}