// Nothing to change here
#pragma once

class NoisyPauli {
   private:
   public:
    static void SingleQubitDepolarise_r(std::vector<double>& cArr, const int numQubits, const int targetQubit, const double err);        // Single-qubit depolarising channel.
    static void NCX_r(std::vector<double>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const double err);  // Two-qubit depolarising channel.
    static void NCZ_r(std::vector<double>& cArr, const int numQubits, const int controlQubit, const int targetQubit, const double err);  // Two-qubit depolarising channel.
};