#include "QubitCode.cuh"

void QubitCode::applyCorrection_flag(cuda_global& dummy, cuda_global& input, unsigned i, unsigned flag, dim3 gridSize, cudaStream_t stream) {
    for (int j = 0; j < NUMQUBITS; j++)
    {
        switch (corrections_flag[flag * (DIM_ANCILLA - 1) + i][j])
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

void QubitCode::Sn_flag(cuda_global& dummy, cuda_global& input, const unsigned num_anc, const unsigned numQubits, const datatype err, const datatype errTeleport, const datatype ancillaFlagErr,
                        std::vector<cudaStream_t>& streams, std::vector<cuda_global>& d_temps) {
    const unsigned N = 1 << numQubits;
    dim3 gridSize = getGridSize(N, N);

    unsigned j = 0;
    for (unsigned i = 2; i < numQubits; i++) {
        if (j == 1 || j == STABILIZER_GATES) {
            NoisyPauli::NCG_dm('X', dummy, input, N, 0, 1, ancillaFlagErr, streams, d_temps);
        } else {
            NoisyPauli::NCG_dm(stabilizers[num_anc][i - 2], dummy, input, N, 0, i, (j == 2) ? err : errTeleport, streams, d_temps);
        }
        swap(input, dummy);
        j++;
    }
}

void QubitCode::OneRoundSequential_flag(cuda_global& input, std::vector<std::vector<cuda_global>>& ORS_DMs_sets, std::vector<cuda_global>& DMs_Nplus2, std::vector<cuda_global>& DMs_Nplus1, std::vector<std::vector<cuda_global>>& d_temps, std::vector<std::vector<cudaStream_t>>& streams_sets, std::vector<std::vector<cublasHandle_t>>& handles_sets, 
                                        const datatype err, const datatype errTeleport, const datatype ancillaFlagErr) {
    auto gridSize_Nplus2 = getGridSize(DIM * 4, DIM * 4);
    auto gridSize_N = getGridSize(DIM, DIM);

    // datatype* DM_H2P0_ptr = nullptr;
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H2P0_ptr, DM_H2P0));

    kron_product<<<gridSize_Nplus2, blockSize, 0, streams_sets[0][0]>>>(DM_H2P0.get(), 4, 4, input.get(), DIM, DIM, DMs_Nplus2[0].get());

    // SET input to 0 HERE
    input.SetZero();

    for (unsigned num_anc = 1; num_anc <= NUM_ANCILLA; num_anc++) {
        QubitCode::Sn_flag(DMs_Nplus2[1], DMs_Nplus2[0], num_anc, NUMQUBITS + 2, err, errTeleport, ancillaFlagErr, streams_sets[0], d_temps[0]);

        CHECK_CUDA(cudaDeviceSynchronize());

        QubitCode::applyAncilla(ORS_DMs_sets[0][1], ancilla_flag[1], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][1]);

        QubitCode::applyAncilla(ORS_DMs_sets[1][2], ancilla_flag[2], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][2]);
        QubitCode::applyAncilla(ORS_DMs_sets[1][3], ancilla_flag[3], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][3]);

        CHECK_CUDA(cudaDeviceSynchronize());

        CHECK_CUBLAS(cublasSaxpy(handles_sets[1][2], DM_SIZE, &alpha, ORS_DMs_sets[1][3].get(), 1, ORS_DMs_sets[1][2].get(), 1));

        CHECK_CUDA(cudaDeviceSynchronize());

        std::jthread t1([&] {
            QubitCode::OneRoundSequential(DMs_Nplus1[1], ORS_DMs_sets[0][1], ORS_DMs_sets[0], DMs_Nplus1[0], d_temps[0], streams_sets[0], handles_sets[0][0], err, errTeleport, QubitCode::correctionMode::noLossW1, 0);
        });

        std::jthread t2([&] {
            QubitCode::OneRoundSequential(DMs_Nplus1[3], ORS_DMs_sets[1][2], ORS_DMs_sets[1], DMs_Nplus1[2], d_temps[1], streams_sets[1], handles_sets[1][0], err, errTeleport, QubitCode::correctionMode::noLossW2, num_anc - 1);

        });

        t1.join();
        t2.join();

        CHECK_CUDA(cudaDeviceSynchronize());
        // Adding results into (zeroed) input matrix
        combine_2d<<<gridSize_N, blockSize, 0, streams_sets[0][0]>>>(input.get(), ORS_DMs_sets[0][1].get(), ORS_DMs_sets[1][2].get(), DIM);
        // CHECK_CUDA(cudaDeviceSynchronize()); 
    }

    QubitCode::applyAncilla(ORS_DMs_sets[0][0], ancilla_flag[0], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][1]);

    CHECK_CUBLAS(cublasSaxpy(handles_sets[0][1], DM_SIZE, &alpha, ORS_DMs_sets[0][0].get(), 1, input.get(), 1));
}