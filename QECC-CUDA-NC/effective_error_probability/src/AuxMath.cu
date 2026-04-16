// Nothing to change here
#include "AuxMath.cuh"

const Matrix I_h(2, 2, {1, 0, 0, 1});
const Matrix X_h(2, 2, {0, 1, 1, 0});
const Matrix Z_h(2, 2, {1, 0, 0, -1});
const Matrix Y_h(2, 2, {0, -1, 1, 0});

// Column-major Kronecker product
// A: m n (col-major), B: p q (col-major)
// Result: (m*p)   (n*q) in col-major order
// Kronecker product (column-major)
Matrix kronecker_col_major(const Matrix& A, const Matrix& B) {
    const unsigned m = A.row, n = A.col;
    const unsigned p = B.row, q = B.col;

    assert(A.mat.size() == m * n);
    assert(B.mat.size() == p * q);

    Matrix result(m * p, n * q);
    std::vector<datatype>& out = result.mat;

    for (unsigned jA = 0; jA < n; ++jA) {
        for (unsigned jB = 0; jB < q; ++jB) {
            unsigned jOut = jA * q + jB;
            for (unsigned iA = 0; iA < m; ++iA) {
                for (unsigned iB = 0; iB < p; ++iB) {
                    unsigned iOut = iA * p + iB;
                    out[jOut * result.row + iOut] =
                        A.mat[jA * m + iA] * B.mat[jB * p + iB];
                }
            }
        }
    }

    return result;
}

double AuxMath::nchoosek(const int n, const int k)
{   // Copied from https://stackoverflow.com/questions/55421835/c-binomial-coefficient-is-too-slow
    // nCk - binomial coeffs
    assert(n >= k);
    double tmp = n - k + 1;
    for (int i = 1; i < k; ++i)
    {
        tmp = tmp * (n - k + 1 + i) / (i + 1);
    }
    return tmp;
}

std::vector<int> AuxMath::AP(const int start, const int end, const int d)
{ 
    // Inclusive of 'end'.
    // AP from 'start' to 'end'
    // E.g. Inputs: start = 0, end = 10, d = 2. Output: {0, 2, 4, 6, 8, 10}
    assert(d > 0);
    const int n = (end - start) / d + 1; // equivalent to std::ceil((end-start)/d).
    std::vector<int> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = start + i * d;
    }
    return out;
}


std::vector<double> AuxMath::linspace(const double start, const double end, const int n)
{
    std::vector<double> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = f_(i, start, end, n);
    }
    return out;
}


std::vector<double> AuxMath::logspace(const double start, const double end, const int n)
{
    std::vector<double> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = std::pow(10.0, f_(i, start, end, n));
    }
    return out;
}


template <typename T>
void AuxMath::print1DVec(std::vector<T> const &v)
{
    std::cout << "\033[0;34m";
    for (auto const &el : v)
    {
        std::cout << el << ", ";
    }
    std::cout << "\033[0m\n";
}