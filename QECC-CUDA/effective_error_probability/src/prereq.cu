#include "prereq.cuh"

std::vector<std::string> corrections;
std::vector<std::string> corrections_erasure;
std::vector<std::string> corrections_flag;
std::vector<std::string> stabilizers;

dim3 blockSize(BLOCK_DIM, BLOCK_DIM);

// __constant__ datatype DM_H1P[4], DM_H1N[4], DM_H2P0[16];
cuda_global DM_H1P(4), DM_H2P0(16);
// std::vector<datatype*> DM_H1(2);
// datatype* DM_H2P0_ptr = nullptr;

std::vector<cuda_global> ancilla, ancilla_flag;

cuda_global::cuda_global(size_t size) {
    CHECK_CUDA(cudaMalloc(&arr, size * sizeof(datatype)));
    this->size = size;
}

// cuda_global::cuda_global(size_t size, datatype *host_arr) : cuda_global(size) {
//     CHECK_CUDA(cudaMemcpy(arr, host_arr, size * sizeof(datatype), cudaMemcpyHostToDevice)); // Note that this is synchronous only
// }

cuda_global::cuda_global(const std::vector<datatype>& host_vec) : cuda_global(host_vec.size()) {
    CHECK_CUDA(cudaMemcpy(arr, host_vec.data(), size * sizeof(datatype), cudaMemcpyHostToDevice));
}

// cuda_global::cuda_global(Dynamic_t &host_mat) : cuda_global(host_mat.size()) {
//     CHECK_CUDA(cudaMemcpy(arr, host_mat.data(), size * sizeof(datatype), cudaMemcpyHostToDevice));
// }

// cuda_global::cuda_global(const cuda_global &other) : cuda_global(other.size) {
//     CHECK_CUDA(cudaMemcpy(arr, other.arr, size * sizeof(datatype), cudaMemcpyDeviceToDevice));
// }

cuda_global::cuda_global(cuda_global &&other) : arr(other.arr), size(other.size) {
    other.arr = nullptr;  // Prevent double free
    other.size = 0;
}

cuda_global &cuda_global::operator=(const cuda_global &other) {
    if (this != &other) {
        if (size != other.size) {
            CHECK_CUDA(cudaFree(arr));
            size = other.size;
            CHECK_CUDA(cudaMalloc(&arr, size * sizeof(datatype)));
        }
        CHECK_CUDA(cudaMemcpy(arr, other.arr, size * sizeof(datatype), cudaMemcpyDeviceToDevice));
    }
    return *this;
}

cuda_global &cuda_global::operator=(cuda_global &&other) {
    if (this != &other) {
        if (arr) {
            CHECK_CUDA(cudaFree(arr));
        }
        arr = other.arr;
        size = other.size;
        other.arr = nullptr;  // Prevent double free
        other.size = 0;
    }
    return *this;
}

cuda_global::~cuda_global() {
    // std::cout << "Destructor called\n";
    if (arr) {
        CHECK_CUDA(cudaFree(arr));
        arr = nullptr;
    }
}

datatype *cuda_global::get(void) const {
    return arr;
}

void cuda_global::SetZero(void) {
    CHECK_CUDA(cudaMemset(arr, 0, size * sizeof(datatype)));
}


std::ostream &operator<<(std::ostream &os, const cuda_global &d_arr) {
    std::vector<datatype> h_arr(d_arr.size);  // Create a host array to store the copied data
    CHECK_CUDA(cudaMemcpy(h_arr.data(), d_arr.arr, d_arr.size * sizeof(datatype), cudaMemcpyDeviceToHost));
    os << "[ ";
    for (auto val : h_arr) {  
        os << val << ' ';
    }
    os << ']';

    return os;
}

void swap(cuda_global &first, cuda_global &second) {
    using std::swap;
    swap(first.arr, second.arr);
    swap(first.size, second.size);
}

void cuda_global::copy_to_host(std::vector<datatype>& host_vec) const {
    if (host_vec.size() != size) {
        host_vec.resize(size);
    }
    CHECK_CUDA(cudaMemcpy(host_vec.data(), arr, size * sizeof(datatype), cudaMemcpyDeviceToHost));
}

void SyncStreams(std::vector<cudaStream_t>& streams) {
    for (auto &stream : streams) {
        CHECK_CUDA(cudaStreamSynchronize(stream));
    }

    // std::for_each(streams.begin(), streams.end(), [](auto &stream){CHECK_CUDA(cudaStreamSynchronize(stream));});
}

dim3 getGridSize(const unsigned mat_rows, const unsigned mat_cols) {
    return dim3((mat_cols + BLOCK_DIM - 1) / BLOCK_DIM, (mat_rows + BLOCK_DIM - 1) / BLOCK_DIM);
}

void readLinesFromFile(std::vector<std::string>& lines, const std::string &filename)    // Modified to not make a large copy of lines when returning
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::copy(std::istream_iterator<std::string>(file), std::istream_iterator<std::string>(), std::back_inserter(lines));
}