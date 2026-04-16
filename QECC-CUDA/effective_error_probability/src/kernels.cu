#include "kernels.cuh"

__global__ void kron_product(datatype* A, int m, int n,
                             datatype* B, int p, int q,
                             datatype* C) {
    int row = blockIdx.y * blockDim.y + threadIdx.y;  // [0, m*p)           // Here row increments before col... are these formulae correct?
    int col = blockIdx.x * blockDim.x + threadIdx.x;  // [0, n*q)

    if (row < m * p && col < n * q) {
        C[row * (n * q) + col] = A[(row / p) * n + (col / q)] * B[(row % p) * q + (col % q)];
    }
}


// __global__ void combine_and_scale_2d(datatype *C, datatype *X, datatype *Y, datatype *Z,
//                                      datatype alpha, datatype beta, unsigned width) {
//     unsigned row = blockIdx.y * blockDim.y + threadIdx.y;
//     unsigned col = blockIdx.x * blockDim.x + threadIdx.x;
//     unsigned i = row * width + col;
//     if (row < width && col < width) {
//         datatype sum = X[i] + Y[i] + Z[i];
//         C[i] = alpha * C[i] + beta * sum;
//     }
// }

__global__ void combine_2d(datatype *C, datatype *X, datatype *Y, unsigned width) {
    unsigned row = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned col = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned i = row * width + col;
    if (row < width && col < width) {
        C[i] += X[i] + Y[i];
    }
}

__global__ void combine_and_scale(datatype *C, datatype **inputs, unsigned num_inputs, 
                                  datatype alpha, datatype beta, unsigned width) {
    unsigned row = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned col = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned i = row * width + col;
    
    if (row < width && col < width) {
        datatype sum = 0;
        for (int j = 0; j < num_inputs; ++j) {
            sum += inputs[j][i];
        }
        C[i] = alpha * C[i] + beta * sum;
    }
}

__global__ void combine_and_scale_3(datatype *C, datatype *X, datatype *Y, datatype *Z, 
                                    datatype alpha, datatype beta, unsigned width) {
    unsigned row = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned col = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned i = row * width + col;
    
    if (row < width && col < width) {
        C[i] = alpha * C[i] + beta * (X[i] + Y[i] + Z[i]);
    }
}

__global__ void combine_and_scale_12(datatype *C, datatype *d1, datatype *d2, datatype *d3, datatype *d4, datatype *d5, datatype *d6, datatype *d7, datatype *d8, datatype *d9, datatype *d10, datatype *d11, datatype *d12, 
                                    datatype alpha, datatype beta, unsigned width) {
    unsigned row = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned col = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned i = row * width + col;
    
    if (row < width && col < width) {
        C[i] = alpha * C[i] + beta * (d1[i] + d2[i] + d3[i] + d4[i] + d5[i] + d6[i] + d7[i] + d8[i] + d9[i] + d10[i] + d11[i] + d12[i]);
    }
}

__global__ void X_dm(datatype* output, datatype* input, unsigned k, unsigned N) {
    // Calculate the mask for flipping the k-th bit
    unsigned mask = (1U << k);
    
    // Calculate global indices: 
    //   row: from x dimension (blockIdx.x and threadIdx.x)
    //   col: from y dimension (blockIdx.y and threadIdx.y)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;

    // Check bounds
    if (row < N && col < N) {
        // Flip the k-th bit of row and column indices
        unsigned flipped_row = row ^ mask;
        unsigned flipped_col = col ^ mask;
        
        // Column-major: index = row + col * N
        output[row + col * N] = input[flipped_row + flipped_col * N];
    }
}

__global__ void Y_dm(datatype* output, datatype* input, unsigned k, unsigned N) {
    // Calculate the mask for flipping the k-th bit
    unsigned mask = (1U << k);
    
    // Calculate global indices (column-major order)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;

    // Check bounds
    if (row < N && col < N) {
        // Flip the k-th bit of row and column indices
        unsigned flipped_row = row ^ mask;
        unsigned flipped_col = col ^ mask;
        unsigned idx = flipped_row + flipped_col * N;

        // Extract k-th bit from row and column indices
        unsigned row_bit = (flipped_row >> k) & 1;
        unsigned col_bit = (flipped_col >> k) & 1;
        
        // Apply phase flip based on the k-th bit
        output[row + col * N] = (row_bit ^ col_bit) ? -input[idx] : input[idx];
    }
}


__global__ void Z_dm(datatype* output, datatype* input, unsigned k, unsigned N) {
    // Calculate global indices (column-major order)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (row < N && col < N) {
        // Calculate column-major index
        unsigned idx = row + col * N;
        
        // Extract k-th bit from row and column indices
        unsigned row_bit = (row >> k) & 1;
        unsigned col_bit = (col >> k) & 1;
        
        // Apply phase flip
        output[idx] = (row_bit ^ col_bit) ? -input[idx] : input[idx];
    }
}

__global__ void CX_dm(datatype* output, datatype* input, 
                      unsigned control_qubit, unsigned target_qubit, unsigned N) {
    // Calculate global indices (column-major order)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (row < N && col < N) {
        // Extract control and target bits
        bool row_control = (row >> control_qubit) & 1;
        bool col_control = (col >> control_qubit) & 1;
        
        // Calculate source indices
        unsigned src_row = row;
        unsigned src_col = col;
        
        // Flip target bit if control is set
        if (row_control) src_row ^= (1U << target_qubit);
        if (col_control) src_col ^= (1U << target_qubit);
        
        // Calculate column-major indices     
        output[row + col * N] = input[src_row + src_col * N];
    }
}

__global__ void CY_dm(datatype* output, datatype* input, 
                      unsigned control_qubit, unsigned target_qubit, unsigned N) {
    // Calculate global indices (column-major order)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (row < N && col < N) {
        // Extract control and target bits from row and column indices
        unsigned control_mask = (1U << control_qubit);
        unsigned target_mask = (1U << target_qubit);
        
        bool row_control = (row & control_mask) != 0;
        bool col_control = (col & control_mask) != 0;
        bool row_target = (row & target_mask) != 0;
        bool col_target = (col & target_mask) != 0;

        // Calculate source indices
        unsigned src_row = row;
        unsigned src_col = col;
        
        // Flip target bit if control is set (for both row and column)
        if (row_control) src_row ^= target_mask;
        if (col_control) src_col ^= target_mask;
        
        // Calculate sign changes:
        // 1. For row: -1 if control set and original target was 0
        // 2. For column: -1 if control set and original target was 1

        // Calculate column-major indices
        unsigned src_idx = src_row + src_col * N;
        unsigned dst_idx = row + col * N;
        
        output[dst_idx] = ~(row_target ^ col_target) ? -input[src_idx] : input[src_idx];

        // float sign = 1.0f;
        // if (row_control) {
        //     if (!row_target) sign *= -1.0f;  // Original target was 0 before flip
        // }
        // if (col_control) {
        //     if (col_target) sign *= -1.0f;   // Original target was 1 before flip
        // }
        // output[dst_idx] = sign * input[src_idx];
    }
}

__global__ void CZ_dm(datatype* output, datatype* input, 
                      unsigned control_qubit, unsigned target_qubit, unsigned N) {
    // Calculate global indices (column-major order)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (row < N && col < N) {
        // Calculate column-major index
        unsigned idx = row + col * N;
        
        // Extract control and target bits
        bool row_control = (row >> control_qubit) & 1;
        bool row_target = (row >> target_qubit) & 1;
        bool col_control = (col >> control_qubit) & 1;
        bool col_target = (col >> target_qubit) & 1;
        
        // Calculate sign based on both control and target bits
        bool row_sign = row_control && row_target;
        bool col_sign = col_control && col_target;
        
        // Apply phase flip
        output[idx] = (row_sign ^ col_sign) ? -input[idx] : input[idx];
    }
}

__global__ void trace_dm(datatype* output, datatype* input, unsigned N) {
    unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < N) {
        atomicAdd(output, input[idx * (idx + 2)]);
    }
}

__global__ void fidel_dm(datatype* output, datatype* input, datatype* logDM, unsigned N) {
    // Calculate global indices (column-major order)
    unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (row < N && col < N) {
        // Calculate column-major index
        unsigned idx = row + col * N;
        atomicAdd(output, input[idx] * logDM[idx]);
    }
}

void NoisyPauli::NCG_dm(char G, cuda_global& output, cuda_global& input, const unsigned N, const unsigned controlQubit, const unsigned targetQubit, const datatype err,
                        std::vector<cudaStream_t>& streams, std::vector<cuda_global>& d_temps) {

    dim3 gridSize = getGridSize(N, N);

    switch (G) {
        case 'I':
            return;
        case 'X':
            CX_dm<<<gridSize, blockSize, 0, streams[0]>>>(output.get(), input.get(), controlQubit, targetQubit, N);
            break;
        case 'Y':
            CY_dm<<<gridSize, blockSize, 0, streams[0]>>>(output.get(), input.get(), controlQubit, targetQubit, N);
            break;
        case 'Z':
            CZ_dm<<<gridSize, blockSize, 0, streams[0]>>>(output.get(), input.get(), controlQubit, targetQubit, N);
            break;
    }

    CUDA_KERNEL_CHECK();

    SyncStreams(streams);

    X_dm<<<gridSize, blockSize, 0, streams[1]>>>(d_temps[0].get(), output.get(), targetQubit, N);
    Y_dm<<<gridSize, blockSize, 0, streams[2]>>>(d_temps[4 + 0].get(), output.get(), targetQubit, N);
    Z_dm<<<gridSize, blockSize, 0, streams[3]>>>(d_temps[8 + 0].get(), output.get(), targetQubit, N);

    CUDA_KERNEL_CHECK();

    SyncStreams(streams);


    X_dm<<<gridSize, blockSize, 0, streams[1]>>>(d_temps[1].get(), d_temps[0].get(), controlQubit, N);
    Y_dm<<<gridSize, blockSize, 0, streams[2]>>>(d_temps[2].get(), d_temps[0].get(), controlQubit, N);
    Z_dm<<<gridSize, blockSize, 0, streams[3]>>>(d_temps[3].get(), d_temps[0].get(), controlQubit, N);

    X_dm<<<gridSize, blockSize, 0, streams[1]>>>(d_temps[4 + 1].get(), d_temps[4].get(), controlQubit, N);
    Y_dm<<<gridSize, blockSize, 0, streams[2]>>>(d_temps[4 + 2].get(), d_temps[4].get(), controlQubit, N);
    Z_dm<<<gridSize, blockSize, 0, streams[3]>>>(d_temps[4 + 3].get(), d_temps[4].get(), controlQubit, N);

    X_dm<<<gridSize, blockSize, 0, streams[1]>>>(d_temps[8 + 1].get(), d_temps[8].get(), controlQubit, N);
    Y_dm<<<gridSize, blockSize, 0, streams[2]>>>(d_temps[8 + 2].get(), d_temps[8].get(), controlQubit, N);
    Z_dm<<<gridSize, blockSize, 0, streams[3]>>>(d_temps[8 + 3].get(), d_temps[8].get(), controlQubit, N);

    CUDA_KERNEL_CHECK();
    
    SyncStreams(streams);

    combine_and_scale_12<<<gridSize, blockSize, 0, streams[0]>>>(output.get(), d_temps[0].get(), d_temps[1].get(), d_temps[2].get(), d_temps[3].get(), d_temps[4].get(), d_temps[5].get(), d_temps[6].get(), d_temps[7].get(), d_temps[8].get(), d_temps[9].get(), d_temps[10].get(), d_temps[11].get(), (1.0f - err) + 3 * err / 15, err / 15, N);
    
    CUDA_KERNEL_CHECK();

    SyncStreams(streams);
}
