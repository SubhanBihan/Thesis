// Nothing to change here
#ifndef _AUXMATHHPP_
#define _AUXMATHHPP_
#include <cstdint>
#include <vector>
#include <cmath>
#include <iostream>
#include <cassert>

#include "params.hpp"

namespace AuxMath{
    dtype nchoosek(const int n, const int k);
    uint64_t nchoosek(uint64_t n, uint64_t k);
    std::vector<int> arange(const int start, const int end, const int d);
    std::vector<dtype> linspace(const dtype start, const dtype end, const int n);
    std::vector<dtype> logspace(const dtype start, const dtype end, const int n);
    dtype f_(const dtype i, const dtype start, const dtype end, const int n);
    template <typename T>
    void print1DVec(std::vector<T> const &v);
}

#endif