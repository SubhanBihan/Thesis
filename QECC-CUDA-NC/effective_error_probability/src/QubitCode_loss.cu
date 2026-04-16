#include "QubitCode.cuh"

void QubitCode::applyCorrection_1erasure(datatype* dummy, datatype* input, unsigned i, unsigned qblost, dim3 gridSize, cudaStream_t stream) {
    for (int j = 0; j < NUMQUBITS; j++)
    {
        switch (corrections_erasure[qblost * (DIM_ANCILLA - 1) + i][j])
        {
        case 'I':
            continue;
        case 'X':
            X_dm<<<gridSize, blockSize, 0, stream>>>(dummy, input, j, NUMQUBITS);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            std::swap(input, dummy);
            break;
        case 'Y':
            Y_dm<<<gridSize, blockSize, 0, stream>>>(dummy, input, j, NUMQUBITS);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            std::swap(input, dummy);
            break;
        case 'Z':
            Z_dm<<<gridSize, blockSize, 0, stream>>>(dummy, input, j, NUMQUBITS);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            std::swap(input, dummy);
            break;
        default:
            break;
        }
    }
}