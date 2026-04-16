#include <iostream>
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>  // For Kronecker/tensor product in Eigen
#include <qpp/qpp.hpp>  // Uncomment if you have qpp installed

// If you have qpp installed:
// #include "qpp.h"

using namespace Eigen;
using namespace std;

int main() {
    cout << "=== EIGEN TENSOR PRODUCT EXAMPLES ===\n\n";
    
    // Example 1: Kronecker Product with Eigen
    cout << "1. Kronecker Product (Tensor Product) with Eigen:\n";
    
    Matrix2d A;  // 2x2 matrix
    A << 1, 2,
         3, 4;
    
    Matrix2d B;  // 2x2 matrix  
    B << 5, 6,
         7, 8;
    
    // Compute Kronecker product A ⊗ B
    auto kron_result = kroneckerProduct(A, B);
    
    cout << "Matrix A:\n" << A << "\n\n";
    cout << "Matrix B:\n" << B << "\n\n";
    cout << "Kronecker Product A ⊗ B:\n" << kron_result << "\n\n";
    
    // Example 2: Tensor product with vectors
    cout << "2. Tensor Product with Vectors:\n";
    
    Vector2d v1(1, 2);
    Vector2d v2(3, 4);
    
    auto vec_kron = kroneckerProduct(v1, v2);
    cout << "Vector v1: " << v1.transpose() << "\n";
    cout << "Vector v2: " << v2.transpose() << "\n";
    cout << "Tensor product v1 ⊗ v2: " << vec_kron.transpose() << "\n\n";
    
    // Example 3: Mixed dimensions
    cout << "3. Mixed Dimensions Example:\n";
    
    Matrix<double, 2, 3> C;  // 2x3 matrix
    C << 1, 2, 3,
         4, 5, 6;
    
    Matrix<double, 2, 2> D;  // 2x2 matrix
    D << 7, 8,
         9, 10;
    
    auto mixed_kron = kroneckerProduct(C, D);
    cout << "Matrix C (2x3):\n" << C << "\n\n";
    cout << "Matrix D (2x2):\n" << D << "\n\n";
    cout << "Kronecker Product C ⊗ D (4x6):\n" << mixed_kron << "\n\n";
    
    cout << "=== QPP LIBRARY EXAMPLES (Quantum Computing) ===\n\n";
    

    // Uncomment if you have qpp installed
    
    using namespace qpp;
    
    cout << "4. qpp Tensor Product Examples:\n";
    
    // qpp uses Eigen matrices internally but provides quantum-specific operations
    
    // Create quantum states (column vectors)
    ket psi0 = st.z0;  // |0⟩ state
    ket psi1 = st.z1;  // |1⟩ state
    
    // Tensor product of quantum states
    ket psi_tensor = kron(psi0, psi1);  // |0⟩ ⊗ |1⟩ = |01⟩
    
    cout << "State |0⟩:\n" << disp(psi0) << "\n";
    cout << "State |1⟩:\n" << disp(psi1) << "\n";
    cout << "Tensor product |0⟩ ⊗ |1⟩:\n" << disp(psi_tensor) << "\n";
    
    // Tensor product of quantum gates
    cmat H = gt.H;     // Hadamard gate
    cmat X = gt.X;     // Pauli-X gate
    
    cmat H_tensor_X = kron(H, X);  // H ⊗ X
    
    cout << "Hadamard gate:\n" << disp(H) << "\n";
    cout << "Pauli-X gate:\n" << disp(X) << "\n";  
    cout << "Tensor product H ⊗ X:\n" << disp(H_tensor_X) << "\n";
    
    // Multiple tensor products
    cmat three_qubit_gate = kron(H, kron(X, gt.Z));  // H ⊗ X ⊗ Z
    cout << "Three-qubit gate H ⊗ X ⊗ Z:\n" << disp(three_qubit_gate) << "\n";

    
    cout << "Note: qpp examples are commented out. To use them:\n";
    cout << "1. Install qpp library from: https://github.com/softwareQinc/qpp\n";
    cout << "2. Uncomment the qpp section above\n";
    cout << "3. Link against qpp in your build system\n";
    
    return 0;
}

/*
Compilation instructions:

For Eigen only:
g++ -std=c++17 -I/path/to/eigen tensor_example.cpp -o tensor_example

For qpp (after installation):
g++ -std=c++17 -I/path/to/qpp/include -I/path/to/eigen tensor_example.cpp -fopenmp -o tensor_example

Key differences:
- Eigen: Uses kroneckerProduct() function, general-purpose linear algebra
- qpp: Uses kron() function, specialized for quantum computing with built-in quantum states and gates
- qpp provides more intuitive quantum mechanics notation and operations
- Both libraries support arbitrary dimensional tensor products
*/