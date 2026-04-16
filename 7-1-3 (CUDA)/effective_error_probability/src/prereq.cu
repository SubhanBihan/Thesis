#include "prereq.cuh"

std::vector<std::string> corrections;
std::vector<std::string> corrections_erasure;
std::vector<std::string> corrections_flag;

cuda_global::cuda_global(size_t size) {
    datatype *temp_arr = nullptr;
    CHECK_CUDA(cudaMalloc(&temp_arr, size * sizeof(datatype)));  // Allows async memcpy later
    arr.reset(temp_arr, [](datatype *ptr) { cudaFree(ptr); });  // Use a custom deleter to free the device memory
    
    this->size = size;
}

cuda_global::cuda_global(size_t size, datatype *host_arr) : cuda_global(size) {
    CHECK_CUDA(cudaMemcpy(arr.get(), host_arr, size * sizeof(datatype), cudaMemcpyHostToDevice)); // Note that this is synchronous only
}

cuda_global::cuda_global(Eigen::Matrix<datatype, Eigen::Dynamic, Eigen::Dynamic> &host_arr) : cuda_global(host_arr.size()) {
    CHECK_CUDA(cudaMemcpy(arr.get(), host_arr.data(), size * sizeof(datatype), cudaMemcpyHostToDevice));  // BEWARE: Column-major format
}

// cuda_global::cuda_global(const cuda_global &other) : cuda_global(other.size) {
//     CHECK_CUDA(cudaMemcpy(arr, other.arr, size * sizeof(datatype), cudaMemcpyDeviceToDevice));
// }

std::ostream &operator<<(std::ostream &os, const cuda_global &d_arr) {
    std::vector<datatype> h_arr(d_arr.size);  // Create a host array to store the copied data
    CHECK_CUDA(cudaMemcpy(h_arr.data(), d_arr.arr.get(), sizeof(datatype), cudaMemcpyDeviceToHost));
    os << "[ ";
    for (auto val : h_arr) {  
        os << val << ' ';
    }
    os << ']';

    return os;
}

void copyUpperTriangularColumnMajor(const Eigen::Matrix<datatype, dynamic, dynamic>& mat, datatype* dest, unsigned reduced_size) {
    int index = 0;
    int n = mat.rows();

    datatype *tmp = (datatype *)malloc(reduced_size * sizeof(datatype));

    auto view = mat.triangularView<Eigen::Upper>();
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i <= j; ++i) {
            tmp[index++] = view(i, j);
        }
    }

    CHECK_CUDA(cudaMemcpy(dest, tmp, reduced_size * sizeof(datatype), cudaMemcpyHostToDevice));
    free(tmp);
}