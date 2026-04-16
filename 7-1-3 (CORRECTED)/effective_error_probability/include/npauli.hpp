// Nothing to change here
#ifndef _NPAULIHPP_
#define _NPAULIHPP_

#include "params.hpp"

class NoisyPauli {
   private:
   public:
    static void SingleQubitDepolarise_r(std::vector<dtype>& cArr, const int numQubits, const int targetQubit, const dtype err);        // Single-qubit depolarising channel.
    static void TwoQubitDepolarise_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const dtype err);  // Two-qubit depolarising channel.
    static void NCX_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const dtype err);  // Two-qubit depolarising channel.
    static void NCZ_r(std::vector<dtype>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const dtype err);  // Two-qubit depolarising channel.
};

#endif