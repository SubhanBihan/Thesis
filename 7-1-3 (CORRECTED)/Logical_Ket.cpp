#include <array>
#include <cmath>  // for sqrt
#include <iostream>  // for printing

// Define the normalization factor
const double factor = 1.0 / std::sqrt(8.0);

// Logical |0_L> codeword for the [[7, 1, 3]] Steane code
std::array<double, 128> logical_ket = {0.0};



// Function to print array as comma-separated values in curly braces
void printArray(const std::array<double, 128>& arr) {
    std::cout << "{";
    for (size_t i = 0; i < arr.size(); ++i) {
        std::cout << arr[i];
        if (i < arr.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "}" << std::endl;
}


int main() {
    // Assign coefficients to |0_L> (for specific basis states)
    logical_ket[0] = 1;        // |0000000>
    logical_ket[85] = 1;       // |1010101>
    logical_ket[51] = 1;       // |0110011>
    logical_ket[102] = 1;      // |1100110>
    logical_ket[15] = 1;       // |0001111>
    logical_ket[90] = 1;       // |1011010>
    logical_ket[60] = 1;       // |0111100>
    logical_ket[105] = 1;      // |1101001>

    // Assign coefficients to |1_L> (for specific basis states)
    logical_ket[127] = 2;      // |1111111>
    logical_ket[42] = 2;       // |0101010>
    logical_ket[76] = 2;       // |1001100>
    logical_ket[25] = 2;       // |0011001>
    logical_ket[112] = 2;      // |1110000>
    logical_ket[37] = 2;       // |0100101>
    logical_ket[67] = 2;       // |1000011>
    logical_ket[22] = 2;       // |0010110>
    
    std::cout << "Logical ket:" << std::endl;
    printArray(logical_ket);

    return 0;
}