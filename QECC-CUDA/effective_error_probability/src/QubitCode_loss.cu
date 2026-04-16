#include "QubitCode.cuh"

void QubitCode::applyCorrection_1erasure(cuda_global& dummy, cuda_global& input, unsigned i, unsigned qblost, dim3 gridSize, cudaStream_t stream) {
    for (int j = 0; j < NUMQUBITS; j++)
    {
        switch (corrections_erasure[qblost * (DIM_ANCILLA - 1) + i][j])
        {
        case 'I':
            continue;
        case 'X':
            X_dm<<<gridSize, blockSize, 0, stream>>>(dummy.get(), input.get(), j, NUMQUBITS);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            swap(input, dummy);
            break;
        case 'Y':
            Y_dm<<<gridSize, blockSize, 0, stream>>>(dummy.get(), input.get(), j, NUMQUBITS);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            swap(input, dummy);
            break;
        case 'Z':
            Z_dm<<<gridSize, blockSize, 0, stream>>>(dummy.get(), input.get(), j, NUMQUBITS);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            swap(input, dummy);
            break;
        // default:
        //     break;
        }
    }
}