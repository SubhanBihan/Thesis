// Nothing to change here
#pragma once

#include <cstdint>
#include <vector>

#define f_(i, start, end, n) (start + i * (end - start) / (static_cast<double>(n) - 1))

namespace AuxMath {
    double nchoosek(const int n, const int k);

    std::vector<int> AP_host(const int start, const int end, const int d);
    __global__ void AP_dev(const int start, const int n, const int d, int* out);

    std::vector<double> linspace_host(const double start, const double end, const int n);
    __global__ void linspace_dev(const double start, const double end, const int n, double* out);

    std::vector<double> logspace_host(const double start, const double end, const int n);
    __global__ void logspace_dev(const double start, const double end, const int n, double* out);

    template <typename T>
    void print1DVec(std::vector<T> const &v);
}