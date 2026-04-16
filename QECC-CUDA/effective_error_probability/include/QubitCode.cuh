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

    void Sn_seq(cuda_global& dummy, cuda_global& input, const unsigned num_anc, const unsigned numQubits, const datatype err, const datatype errTeleport,
                std::vector<cudaStream_t>& streams, std::vector<cuda_global>& d_temps);

    void Sn_flag(cuda_global& dummy, cuda_global& input, const unsigned num_anc, const unsigned numQubits, const datatype err, const datatype errTeleport, const datatype ancillaFlagErr,
                 std::vector<cudaStream_t>& streams, std::vector<cuda_global>& d_temps);

    void applyAncilla(cuda_global& output, cuda_global& ancilla_mat, cuda_global& input, unsigned m, unsigned n, cublasHandle_t handle);

    void applyCorrection(cuda_global& dummy, cuda_global& input, unsigned i, dim3 gridSize, cudaStream_t stream);

    void applyCorrection_1erasure(cuda_global& dummy, cuda_global& input, unsigned i, unsigned qblost, dim3 gridSize, cudaStream_t stream);
    
    void applyCorrection_flag(cuda_global& dummy, cuda_global& input, unsigned i, unsigned flag, dim3 gridSize, cudaStream_t stream);

    void ORS_NestedLoops(cuda_global& dummy, const unsigned num_anc, unsigned &corr_ind, std::vector<cuda_global>& ORS_DMs, cuda_global &DM_Nplus1, std::vector<cuda_global>& d_temps, std::vector<cudaStream_t>& streams, cublasHandle_t handle, 
                          dim3 gridSize_N, dim3 gridSize_Nplus1, const datatype err, const datatype errTeleport, const QubitCode::correctionMode correctionMode, const unsigned qblost_or_flag);
                         
    void OneRoundSequential(cuda_global& dummy, cuda_global& input, std::vector<cuda_global>& ORS_DMs, cuda_global &DM_Nplus1, std::vector<cuda_global>& d_temps, std::vector<cudaStream_t>& streams, cublasHandle_t handle, 
                            const datatype err, const datatype errTeleport, const QubitCode::correctionMode correctionMode, const unsigned qblost_or_flag);

    void OneRoundSequential_flag(cuda_global& input, std::vector<std::vector<cuda_global>>& ORS_DMs_sets, std::vector<cuda_global>& DMs_Nplus2, std::vector<cuda_global>& DMs_Nplus1, std::vector<std::vector<cuda_global>>& d_temps, std::vector<std::vector<cudaStream_t>>& streams_sets, std::vector<std::vector<cublasHandle_t>>& handles_sets, 
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