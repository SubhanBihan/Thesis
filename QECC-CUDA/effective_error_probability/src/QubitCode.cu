#include "QubitCode.cuh"

void QubitCode::Sn_seq(cuda_global& dummy, cuda_global& input, const unsigned num_anc, const unsigned numQubits, const datatype err, const datatype errTeleport,
                       std::vector<cudaStream_t>& streams, std::vector<cuda_global>& d_temps) {
    const unsigned N = 1 << numQubits;
    // dim3 gridSize = getGridSize(N, N);
    
    //unsigned j = 0;
    bool j = (num_anc == 1);
    for (unsigned i = 1; i < numQubits; i++) {
        if (stabilizers[num_anc - 1][i - 1] != 'I') {
            //NoisyPauli::NCG_dm(stabilizers[num_anc][i - 1], dummy, input, N, 0, i, (num_anc == 1 && j == 0) ? err : errTeleport, streams, d_temps);
            NoisyPauli::NCG_dm(stabilizers[num_anc - 1][i - 1], dummy, input, N, 0, i, j ? err : errTeleport, streams, d_temps);    // This is ultimately called by flag version too... is it logically correct for that too???
            swap(input, dummy);
            //j++;
            j = false;
        }
    }
}

void QubitCode::applyAncilla(cuda_global& output, cuda_global& ancilla_mat, cuda_global& input, unsigned m, unsigned n, cublasHandle_t handle){
    cuda_global tmp(m * n);
    tmp.SetZero();

    CHECK_CUBLAS(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                             m, n, n, &alpha,
                             ancilla_mat.get(), m,
                             input.get(), n,
                             &beta,
                             tmp.get(), m));
                
    CHECK_CUBLAS(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T,
                             m, m, n, &alpha,
                             tmp.get(), m,
                             ancilla_mat.get(), m,
                             &beta,
                             output.get(), m));
}

void QubitCode::applyCorrection(cuda_global& dummy, cuda_global& input, unsigned i, dim3 gridSize, cudaStream_t stream) {
    for (int j = 0; j < NUMQUBITS; j++)
    {
        switch (corrections[i][j])
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
            //default:
            //    break;
        }
    }
}

void QubitCode::OneRoundSequential(cuda_global& dummy, cuda_global& input, std::vector<cuda_global>& ORS_DMs, cuda_global &DM_Nplus1, std::vector<cuda_global>& d_temps, std::vector<cudaStream_t>& streams, cublasHandle_t handle, 
                                   const datatype err, const datatype errTeleport, const QubitCode::correctionMode correctionMode, const unsigned qblost_or_flag) {
    ORS_DMs[0] = input;

    auto gridSize_N = getGridSize(DIM, DIM);
    auto gridSize_Nplus1 = getGridSize(DIM * 2, DIM * 2);
    const unsigned num_anc = 1;
    unsigned corr_ind = 0;

    ORS_DMs[NUM_ANCILLA + 1].SetZero();

    ORS_NestedLoops(dummy, num_anc, corr_ind, ORS_DMs, DM_Nplus1, d_temps, streams, handle, gridSize_N, gridSize_Nplus1, err, errTeleport, correctionMode, qblost_or_flag);

    input = ORS_DMs[NUM_ANCILLA + 1];
}


void QubitCode::ORS_NestedLoops(cuda_global& dummy, const unsigned num_anc, unsigned &corr_ind, std::vector<cuda_global>& ORS_DMs, cuda_global &DM_Nplus1, std::vector<cuda_global>& d_temps, std::vector<cudaStream_t>& streams, cublasHandle_t handle, 
                                dim3 gridSize_N, dim3 gridSize_Nplus1, const datatype err, const datatype errTeleport, const QubitCode::correctionMode correctionMode, const unsigned qblost_or_flag) {

    for (int i = 0; i < 2; i++) {
        kron_product<<<gridSize_Nplus1, blockSize, 0, streams[0]>>>(DM_H1P.get(), 2, 2, ORS_DMs[num_anc - 1].get(), DIM, DIM, DM_Nplus1.get());
        CUDA_KERNEL_CHECK();
        QubitCode::Sn_seq(dummy, DM_Nplus1, num_anc, NUMQUBITS + 1, err, errTeleport, streams, d_temps);
        QubitCode::applyAncilla(ORS_DMs[num_anc], ancilla[i], DM_Nplus1, DIM, DIM * 2, handle);

        SyncStreams(streams);
        
        if (num_anc == NUM_ANCILLA)
        {
            switch (correctionMode)
            {
                case QubitCode::correctionMode::noLossW1:

                    if(corr_ind) QubitCode::applyCorrection(dummy, ORS_DMs[NUM_ANCILLA], corr_ind - 1, gridSize_N, streams[0]);

                    CHECK_CUBLAS(cublasSaxpy(handle, DM_SIZE, &alpha, ORS_DMs[NUM_ANCILLA].get(), 1, ORS_DMs[NUM_ANCILLA + 1].get(), 1));

                    break;

                case QubitCode::correctionMode::oneLoss:

                    if(corr_ind) QubitCode::applyCorrection_1erasure(dummy, ORS_DMs[NUM_ANCILLA], corr_ind - 1, qblost_or_flag, gridSize_N, streams[0]);

                    CHECK_CUBLAS(cublasSaxpy(handle, DM_SIZE, &alpha, ORS_DMs[NUM_ANCILLA].get(), 1, ORS_DMs[NUM_ANCILLA + 1].get(), 1));

                    break;

                case QubitCode::correctionMode::noLossW2:

                    if(corr_ind) QubitCode::applyCorrection_flag(dummy, ORS_DMs[NUM_ANCILLA], corr_ind - 1, qblost_or_flag, gridSize_N, streams[0]);

                    CHECK_CUBLAS(cublasSaxpy(handle, DM_SIZE, &alpha, ORS_DMs[NUM_ANCILLA].get(), 1, ORS_DMs[NUM_ANCILLA + 1].get(), 1));

                    break;

                // default:
                //     break;
            }

            SyncStreams(streams);
            corr_ind++;
        }
        else
        {
            ORS_NestedLoops(dummy, num_anc + 1, corr_ind, ORS_DMs, DM_Nplus1, d_temps, streams, handle, gridSize_N, gridSize_Nplus1, err, errTeleport, correctionMode, qblost_or_flag);
        }

        // std::cout << "2 - HERE!!!" << std::endl;
        
    }
    SyncStreams(streams);
    // std::cout << "3 - HERE!!!" << std::endl;
}
