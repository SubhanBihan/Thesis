// Nothing to change here
#include <cstdint>
#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>

#include "AuxMath.cuh"

double AuxMath::nchoosek(const int n, const int k)
{ // Copied from https://stackoverflow.com/questions/55421835/c-binomial-coefficient-is-too-slow
    // nCk - binomial coeffs
    assert(n >= k);
    double tmp = n - k + 1;
    for (int i = 1; i < k; ++i)
    {
        tmp = tmp * (n - k + 1 + i) / (i + 1);
    }
    return tmp;
}

std::vector<int> AuxMath::AP_host(const int start, const int end, const int d)
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

__global__ void AuxMath::AP_dev(const int start, const int n, const int d, int* out)
{   
    // Unfortunately kernel (with mem transfer OH) is slower!
    // So only use if keeping the AP on device
    // Recommended: <<<ceil(n/256), 256>>> where n = (end-start)/d + 1                                                                             
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < n) out[i] = start + i * d; 
}

std::vector<double> AuxMath::linspace_host(const double start, const double end, const int n)
{
    std::vector<double> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = f_(i, start, end, n);
    }
    return out;
}

__global__ void AuxMath::linspace_dev(const double start, const double end, const int n, double* out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < n) out[i] = f_(i, start, end, n);
}

std::vector<double> AuxMath::logspace_host(const double start, const double end, const int n)
{
    std::vector<double> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = std::pow(10.0, f_(i, start, end, n));
    }
    return out;
}

__global__ void AuxMath::logspace_dev(const double start, const double end, const int n, double* out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < n) out[i] = std::pow(10.0, f_(i, start, end, n));  // Doesn't seem like there's any good API alt
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
template void AuxMath::print1DVec(std::vector<uint16_t> const &v);
template void AuxMath::print1DVec(std::vector<int> const &v);
template void AuxMath::print1DVec(std::vector<double> const &v);