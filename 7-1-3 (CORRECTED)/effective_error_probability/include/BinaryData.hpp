// Generates binary files - could consider as 'printing' the result data
// Minor changes related to number of encoded qubits must be made
// Changes already made here

#ifndef _BINARYDATAHPP_
#define _BINARYDATAHPP_

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "params.hpp"

namespace BinaryData
{
    struct ECFlagInfo
    {
        // std::array<std::vector<dtype>, NUMQUBITS> traces; // These act like probabilities.
        std::array<dtype, NUMPOINTS> errs0;
        std::array<dtype, NUMPOINTS> errs1;
        void read(std::string const &pathFile0, std::string const &pathFile1)
        {
            std::ifstream file0(pathFile0, std::ios::binary);

            file0.seekg(0, std::ios::beg);
            // file0.ignore((2 * NUMQUBITS) * 8 * NUMPOINTS); // Reference: https://stackoverflow.com/questions/18640001/read-several-bytes-jump-over-n-bytes-and-then-read-several-bytes-again-how

            file0.read(reinterpret_cast<char *>(&errs0[0]), errs0.size() * sizeof(dtype)); // Reference: https://stackoverflow.com/questions/15143670/how-can-i-use-fread-on-a-binary-file-to-read-the-data-into-a-stdvector
            file0.close();

            std::ifstream file1(pathFile1, std::ios::binary);

            file1.seekg(0, std::ios::beg);
            // file1.ignore(NUMQUBITS * 8 * NUMPOINTS); // Reference: https://stackoverflow.com/questions/18640001/

            // for (int i = 0; i < NUMQUBITS; i++)     // Made a loop.
            // {
            //     traces[i].resize(NUMPOINTS);
            //     file1.read(reinterpret_cast<char *>(&traces[i][0]), traces[i].size() * sizeof(dtype));
            // }
            
            
            file1.read(reinterpret_cast<char *>(&errs1[0]), errs1.size() * sizeof(dtype)); // Reference: https://stackoverflow.com/questions/15143670/how-can-i-use-fread-on-a-binary-file-to-read-the-data-into-a-stdvector
            file1.close();
        }
        void println(std::array<dtype, NUMPOINTS> const &array) const
        {
            for (size_t i = 0; i < NUMPOINTS; i++)
            {
                std::cout << std::setprecision(std::numeric_limits<dtype>::digits10) << array[i] << '\n';
            }
        }

        // Pg. 24-25 of Masters' Thesis
        dtype getRecursiveFidelity(const size_t n, const uint64_t c) const
        {
            dtype denom;
            const dtype fidelity0 = 1 - errs0[n];
            const dtype fidelity1 = 1 - errs1[n];
            const dtype zeta = std::sqrt((NUMQUBITS - 1) * fidelity1 - (NUMQUBITS - 2) * (fidelity0 * fidelity0));
            switch (c)
            {
            case UINT64_C(1):
                return fidelity0;
                break;
            case UINT64_C(2):
                return fidelity1;
                break;
            default:
                denom = std::pow(2.0, c + 1);

                // So fidelity0 = alpha_1
                if (!std::isnan(denom))
                {
                    return (std::pow(fidelity0 + zeta, c + 1) - std::pow(fidelity0 - zeta, c + 1)) / zeta / denom;
                }
                else
                {
                    return 0.0;
                }
                break;
            }
        }
        void clear()
        {
            // for (size_t i = 0; i < NUMQUBITS; i++)
            // {
            //     if (!traces[i].empty())
            //     {
            //         traces[i].clear();
            //     }
            // }
            errs0.fill(0);
            errs1.fill(0);
        }
    };
    
    struct ECErasureInfo
    {
        std::array<dtype, NUMPOINTS> errs;
        void read(std::string const &pathFile)
        {
            std::ifstream input_file;
            if (std::filesystem::exists(pathFile))
            {
                input_file.open(pathFile, std::ios::binary);

                input_file.seekg(0, std::ios::beg);
                // input_file.ignore((2 * NUMQUBITS) * 8 * NUMPOINTS); // Reference: https://stackoverflow.com/questions/18640001/read-several-bytes-jump-over-n-bytes-and-then-read-several-bytes-again-how
                
                input_file.read(reinterpret_cast<char *>(&this->errs[0]), this->errs.size() * sizeof(dtype)); // Reference: https://stackoverflow.com/questions/15143670/how-can-i-use-fread-on-a-binary-file-to-read-the-data-into-a-stdvector
                input_file.close();
            }
            else
            {
                std::cerr << "\x1b[0;91mFile not found!\n";
            }
        }
    };

} // namespace BinaryData

#endif