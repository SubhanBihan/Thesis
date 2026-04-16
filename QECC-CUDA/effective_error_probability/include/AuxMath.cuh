#pragma once

#include <cstdint>
#include <vector>
#include <ranges>
#include <iostream>
#include <cassert>

#include "params.cuh"

#define f_(i, start, end, n) (start + i * (end - start) / (static_cast<double>(n) - 1))

struct Matrix {
    unsigned row{}, col{};
    std::vector<datatype> mat;

    // Size constructor
    Matrix(unsigned r, unsigned c) : row(r), col(c), mat(r * c) {}

    // Constructor with dimensions and initializer list
    Matrix(unsigned r, unsigned c, std::initializer_list<datatype> init)
        : row(r), col(c), mat(init)
    {
        assert(init.size() == r * c && "Initializer list size must match matrix dimensions.");
    }
};

Matrix kronecker_col_major(const Matrix& A, const Matrix& B);

namespace AuxMath {
    double nchoosek(const int n, const int k);

    std::vector<int> AP(const int start, const int end, const int d);

    std::vector<double> linspace(const double start, const double end, const int n);

    std::vector<double> logspace(const double start, const double end, const int n);

    template <typename T>
    void print1DVec(std::vector<T> const &v);
}

extern const Matrix I_h, X_h, Z_h, Y_h;

