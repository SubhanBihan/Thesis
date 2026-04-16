// Nothing to change here
#include <cstdint>
#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>

#include "AuxMath.hpp"

dtype AuxMath::nchoosek(const int n, const int k)
{ // Copied from https://stackoverflow.com/questions/55421835/c-binomial-coefficient-is-too-slow
    // nCk - binomial coeffs
    assert(n >= k);
    dtype tmp = n - k + 1;
    for (int i = 1; i < k; ++i)
    {
        tmp = tmp * (n - k + 1 + i) / (i + 1);
    }
    return tmp;
}

dtype AuxMath::f_(const dtype i, const dtype start, const dtype end, const int n)
{
    return start + i * (end - start) / (static_cast<dtype>(n) - 1);
}

std::vector<int> AuxMath::arange(const int start, const int end, const int d)
{ // Inclusive of 'end'.
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

std::vector<dtype> AuxMath::linspace(const dtype start, const dtype end, const int n)
{
    std::vector<dtype> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = AuxMath::f_(i, start, end, n);
    }
    return out;
}

std::vector<dtype> AuxMath::logspace(const dtype start, const dtype end, const int n)
{
    std::vector<dtype> out(n);
    for (int i = 0; i < n; i++)
    {
        out[i] = std::pow(10.0, AuxMath::f_(i, start, end, n));
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
template void AuxMath::print1DVec(std::vector<uint16_t> const &v);
template void AuxMath::print1DVec(std::vector<int> const &v);
template void AuxMath::print1DVec(std::vector<dtype> const &v);