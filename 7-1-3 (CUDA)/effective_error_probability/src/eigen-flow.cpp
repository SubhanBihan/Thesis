#include <iostream>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include <cassert>

using namespace Eigen;
using namespace std;
using Complex = complex<double>;

class PackedDensityMatrix {
private:
    int numQubits_;
    int dim_;  // 2^numQubits
    vector<Complex> packedData_;
    
    // Convert (i,j) matrix indices to packed storage index
    size_t getPackedIndex(int i, int j) const {
        // Store upper triangle in column-major order
        // For upper triangle: i <= j
        if (i > j) swap(i, j);  // Ensure i <= j
        
        // Column-major packed index for upper triangle
        return j * (j + 1) / 2 + i;
    }
    
public:
    PackedDensityMatrix(int numQubits) 
        : numQubits_(numQubits), dim_(1 << numQubits) {
        
        // Calculate packed size: n(n+1)/2 for upper triangle including diagonal
        size_t packedSize = (static_cast<size_t>(dim_) * (dim_ + 1)) / 2;
        packedData_.resize(packedSize);
        fill(packedData_.begin(), packedData_.end(), Complex(0.0, 0.0));
    }
    
    // Get element (i,j) - handles Hermitian property
    Complex get(int i, int j) const {
        if (i <= j) {
            return packedData_[getPackedIndex(i, j)];
        } else {
            // For Hermitian matrix: rho(i,j) = conj(rho(j,i))
            return conj(packedData_[getPackedIndex(j, i)]);
        }
    }
    
    // Set element (i,j) - automatically maintains Hermitian property
    void set(int i, int j, const Complex& value) {
        if (i <= j) {
            packedData_[getPackedIndex(i, j)] = value;
        } else {
            // Store conjugate in upper triangle
            packedData_[getPackedIndex(j, i)] = conj(value);
        }
    }
    
    // Convert to full Eigen matrix for operations
    MatrixXcd toFullMatrix() const {
        MatrixXcd full(dim_, dim_);
        for (int i = 0; i < dim_; ++i) {
            for (int j = 0; j < dim_; ++j) {
                full(i, j) = get(i, j);
            }
        }
        return full;
    }
    
    // Load from full matrix (extracts upper triangle)
    void fromFullMatrix(const MatrixXcd& full) {
        assert(full.rows() == dim_ && full.cols() == dim_);
        for (int i = 0; i < dim_; ++i) {
            for (int j = i; j < dim_; ++j) {
                packedData_[getPackedIndex(i, j)] = full(i, j);
            }
        }
    }
    
    // Apply unitary gate: rho' = U * rho * U†
    void applyGate(const MatrixXcd& U) {
        MatrixXcd full = toFullMatrix();
        MatrixXcd result = U * full * U.adjoint();
        fromFullMatrix(result);
    }
    
    // Apply gate to specific qubits (with tensor product expansion)
    void applyGateToQubits(const MatrixXcd& gate, const vector<int>& targetQubits) {
        // For simplicity, convert to full matrix, apply gate, convert back
        // In practice, you might want more efficient implementations
        
        MatrixXcd fullGate = expandGateToFullSpace(gate, targetQubits);
        applyGate(fullGate);
    }
    
    // Expand single/multi-qubit gate to full Hilbert space
    MatrixXcd expandGateToFullSpace(const MatrixXcd& gate, const vector<int>& targetQubits) const {
        // This is a simplified implementation
        // For a complete implementation, you'd want to use more efficient tensor product operations
        
        if (targetQubits.size() == 1) {
            // Single qubit gate
            int target = targetQubits[0];
            MatrixXcd fullGate = MatrixXcd::Identity(dim_, dim_);
            
            // Apply gate to target qubit using Kronecker products
            MatrixXcd I2 = MatrixXcd::Identity(2, 2);
            
            // Build: I ⊗ I ⊗ ... ⊗ gate ⊗ ... ⊗ I
            fullGate = MatrixXcd::Identity(1, 1);
            for (int q = 0; q < numQubits_; ++q) {
                if (q == target) {
                    fullGate = kroneckerProduct(fullGate, gate).eval();
                } else {
                    fullGate = kroneckerProduct(fullGate, I2).eval();
                }
            }
            return fullGate;
        }
        
        // For multi-qubit gates, you'd need more complex logic
        // This is just a placeholder
        return gate;
    }
    
    // Get trace (should be 1 for valid density matrix)
    Complex trace() const {
        Complex tr(0.0, 0.0);
        for (int i = 0; i < dim_; ++i) {
            tr += get(i, i);
        }
        return tr;
    }
    
    // Print packed storage info
    void printInfo() const {
        cout << "Packed Density Matrix Info:\n";
        cout << "  Qubits: " << numQubits_ << "\n";
        cout << "  Dimension: " << dim_ << "x" << dim_ << "\n";
        cout << "  Packed size: " << packedData_.size() << " elements\n";
        cout << "  Memory savings: " << (1.0 - 2.0 * packedData_.size() / (dim_ * dim_)) * 100 << "%\n";
        cout << "  Trace: " << trace() << "\n\n";
    }
    
    int getDimension() const { return dim_; }
    int getNumQubits() const { return numQubits_; }
    const vector<Complex>& getPackedData() const { return packedData_; }
};

// Quantum gate definitions
namespace Gates {
    const MatrixXcd I = (MatrixXcd(2,2) << 1, 0, 0, 1).finished();
    const MatrixXcd X = (MatrixXcd(2,2) << 0, 1, 1, 0).finished();
    const MatrixXcd Y = (MatrixXcd(2,2) << 0, Complex(0,-1), Complex(0,1), 0).finished();
    const MatrixXcd Z = (MatrixXcd(2,2) << 1, 0, 0, -1).finished();
    const MatrixXcd H = (MatrixXcd(2,2) << 1, 1, 1, -1).finished() / sqrt(2.0);
    
    MatrixXcd RZ(double theta) {
        MatrixXcd rz(2, 2);
        rz << exp(Complex(0, -theta/2)), 0, 
              0, exp(Complex(0, theta/2));
        return rz;
    }
}

int main() {
    cout << "=== PACKED DENSITY MATRIX OPERATIONS ===\n\n";
    
    // Example with 3-qubit system (smaller for illustration)
    const int numQubits = 3;
    PackedDensityMatrix rho(numQubits);
    
    cout << "1. INITIALIZATION\n";
    cout << "-----------------\n";
    
    // Initialize to |000⟩⟨000| state (computational basis state)
    rho.set(0, 0, Complex(1.0, 0.0));  // Only (0,0) element is 1
    rho.printInfo();
    
    cout << "Initial state |000⟩⟨000|:\n";
    cout << rho.toFullMatrix() << "\n\n";
    
    
    cout << "2. APPLYING SINGLE QUBIT GATES\n";
    cout << "-------------------------------\n";
    
    // Apply Hadamard gate to qubit 0
    cout << "Applying Hadamard to qubit 0...\n";
    rho.applyGateToQubits(Gates::H, {0});
    
    cout << "After H on qubit 0:\n";
    cout << rho.toFullMatrix() << "\n";
    cout << "Trace: " << rho.trace() << "\n\n";
    
    
    cout << "3. APPLYING ROTATION GATES\n";
    cout << "--------------------------\n";
    
    // Apply RZ rotation to qubit 1
    double theta = M_PI / 4;  // 45 degrees
    cout << "Applying RZ(π/4) to qubit 1...\n";
    rho.applyGateToQubits(Gates::RZ(theta), {1});
    
    cout << "After RZ(π/4) on qubit 1:\n";
    cout << rho.toFullMatrix() << "\n";
    cout << "Trace: " << rho.trace() << "\n\n";
    
    
    cout << "4. MEMORY USAGE COMPARISON\n";
    cout << "--------------------------\n";
    
    // Compare memory usage for different system sizes
    vector<int> qubitCounts = {2, 3, 4, 5, 6};
    
    cout << "Memory usage comparison:\n";
    cout << "Qubits | Full Matrix | Packed Storage | Savings\n";
    cout << "-------|-------------|----------------|--------\n";
    
    for (int nq : qubitCounts) {
        int dim = 1 << nq;
        size_t fullSize = dim * dim * sizeof(Complex);
        size_t packedSize = (static_cast<size_t>(dim) * (dim + 1) / 2) * sizeof(Complex);
        double savings = (1.0 - static_cast<double>(packedSize) / fullSize) * 100;
        
        cout << "   " << nq << "   |   " 
             << fullSize / 1024 << " KB    |     " 
             << packedSize / 1024 << " KB     |  " 
             << fixed << setprecision(1) << savings << "%\n";
    }
    cout << "\n";
    
    
    cout << "5. LARGE SYSTEM EXAMPLE (CONCEPTUAL)\n";
    cout << "------------------------------------\n";
    
    // Show what the storage would look like for larger systems
    cout << "For an 11-qubit system:\n";
    int dim11 = 1 << 11;  // 2048x2048
    size_t fullSize11 = static_cast<size_t>(dim11) * dim11 * sizeof(Complex);
    size_t packedSize11 = (static_cast<size_t>(dim11) * (dim11 + 1) / 2) * sizeof(Complex);
    
    cout << "  Dimension: " << dim11 << "x" << dim11 << "\n";
    cout << "  Full matrix: " << fullSize11 / (1024*1024) << " MB\n";
    cout << "  Packed storage: " << packedSize11 / (1024*1024) << " MB\n";
    cout << "  Memory savings: " << (1.0 - static_cast<double>(packedSize11) / fullSize11) * 100 << "%\n\n";
    
    
    cout << "6. EFFICIENT GATE APPLICATION STRATEGY\n";
    cout << "--------------------------------------\n";
    
    cout << "For very large systems, consider:\n";
    cout << "• Implement gates directly on packed storage when possible\n";
    cout << "• Use sparse matrix representations for gates\n";
    cout << "• Apply gates only to relevant subspaces\n";
    cout << "• Use matrix-free methods for large matrix-vector products\n";
    cout << "• Consider approximate methods for very large systems\n\n";
    
    
    cout << "7. VERIFICATION\n";
    cout << "---------------\n";
    
    // Verify Hermitian property
    MatrixXcd full = rho.toFullMatrix();
    MatrixXcd hermitianCheck = full - full.adjoint();
    double hermitianError = hermitianCheck.norm();
    
    cout << "Hermitian property check (should be ~0): " << hermitianError << "\n";
    cout << "Trace (should be 1): " << rho.trace() << "\n";
    
    // Check positive semidefinite (eigenvalues should be non-negative)
    SelfAdjointEigenSolver<MatrixXcd> eigensolver(full);
    VectorXd eigenvals = eigensolver.eigenvalues();
    double minEigenval = eigenvals.minCoeff();
    
    cout << "Minimum eigenvalue (should be >= 0): " << minEigenval << "\n";
    cout << "Valid density matrix: " << (abs(rho.trace() - Complex(1.0)) < 1e-10 && 
                                         hermitianError < 1e-10 && 
                                         minEigenval >= -1e-10 ? "Yes" : "No") << "\n";
    
    return 0;
}

/*
Key Features of This Implementation:

1. MEMORY EFFICIENT STORAGE:
   - Only stores upper triangle + diagonal
   - ~50% memory savings for large matrices
   - Handles Hermitian property automatically

2. EIGEN INTEGRATION:
   - Converts between packed and full format
   - Uses Eigen for matrix operations
   - Maintains numerical precision

3. QUANTUM GATE OPERATIONS:
   - Single qubit gates with tensor product expansion
   - Maintains trace and Hermitian properties
   - Easy to extend for multi-qubit gates

4. SCALABILITY CONSIDERATIONS:
   - Works with your storage formula
   - Can be optimized further for very large systems
   - Provides foundation for more advanced techniques

For 11-qubit systems and beyond, consider:
- Implementing gates directly on packed storage
- Using GPU acceleration
- Sparse matrix techniques
- Matrix-free iterative methods

Compilation:
g++ -std=c++17 -I/path/to/eigen -O3 packed_dm_example.cpp -o packed_dm_example
*/