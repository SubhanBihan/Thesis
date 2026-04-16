#include "QubitCode.cuh"

void QubitCode::applyCorrection_flag(datatype* dummy, datatype* input, unsigned i, unsigned flag, dim3 gridSize, cudaStream_t stream) {
    for (int j = 0; j < NUMQUBITS; j++)
    {
        switch (corrections_flag[flag * (DIM_ANCILLA - 1) + i][j])
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

void QubitCode::Sn_flag(datatype* dummy, datatype* input, unsigned n, const unsigned numQubits, const datatype err, const datatype errTeleport, const datatype ancillaFlagErr,
                        std::vector<cudaStream_t>& streams, datatype** d_temps) {
    const unsigned N = 1 << numQubits;
    dim3 gridSize = getGridSize(N, N);

    unsigned j = 0;
    for (unsigned i = 2; i < numQubits; i++) {
        if (j == 1 || j == STABILIZER_GATES) {
            NoisyPauli::NCG_dm('X', dummy, input, N, 0, 1, ancillaFlagErr, streams, d_temps);
            std::swap(input, dummy);
            j++;
        }
        else {
            NoisyPauli::NCG_dm(stabilizers[n][i - 2], dummy, input, N, 0, i, (j == 2) ? err : errTeleport, streams, d_temps);
            std::swap(input, dummy);
            j++;

            // switch (stabilizers[n][i - 2]) {
            //     case 'I':
            //         continue;
            //     case 'X':
                    // NoisyPauli::NCX_dm(dummy, input, N, 0, i, (j == 2) ? err : errTeleport, streams, d_temps);
                    // std::swap(input, dummy);
                    // j++;
            //         break;
            //     case 'Z':
            //         NoisyPauli::NCZ_dm(dummy, input, N, 0, i, (j == 2) ? err : errTeleport, streams, d_temps);
            //         std::swap(input, dummy);
            //         j++;
            //         break;
            //     default:
            //         break;
            // }
        }
    }
}

void QubitCode::OneRoundSequential_flag(datatype* input, std::vector<std::vector<datatype*>>& ORS_DMs_sets, std::vector<datatype*>& DMs_Nplus2, std::vector<datatype*>& DMs_Nplus1, std::array<datatype**, 2>& d_temps, std::array<std::vector<cudaStream_t>, 2>& streams_sets, std::array<std::vector<cublasHandle_t>, 2>& handles_sets, 
                                        const datatype err, const datatype errTeleport, const datatype ancillaFlagErr) {
    auto gridSize_Nplus2 = getGridSize(DIM * 4, DIM * 4);
    auto gridSize_N = getGridSize(DIM, DIM);

    // datatype* DM_H2P0_ptr = nullptr;
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H2P0_ptr, DM_H2P0));


    kron_product<<<gridSize_Nplus2, blockSize, 0, streams_sets[0][0]>>>(DM_H2P0, 4, 4, input, DIM, DIM, DMs_Nplus2[0]);

    // SET input to 0 HERE
    // input.SetZero();
    CHECK_CUDA(cudaMemset(input, 0, DM_SIZE * sizeof(datatype)));


    for (unsigned num_anc = 1; num_anc <= NUM_ANCILLA; num_anc++) {
        QubitCode::Sn_flag(DMs_Nplus2[1], DMs_Nplus2[0], num_anc, NUMQUBITS + 2, err, errTeleport, ancillaFlagErr, streams_sets[0], d_temps[0]);

        QubitCode::applyAncilla(ORS_DMs_sets[0][1], ancilla_flag[1], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][1]);

        QubitCode::applyAncilla(ORS_DMs_sets[1][2], ancilla_flag[2], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][2]);
        QubitCode::applyAncilla(ORS_DMs_sets[1][3], ancilla_flag[3], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][3]);

        CHECK_CUDA(cudaDeviceSynchronize());

        CHECK_CUBLAS(cublasSaxpy(handles_sets[1][2], DM_SIZE, &alpha, ORS_DMs_sets[1][3], 1, ORS_DMs_sets[1][2], 1));

        CHECK_CUDA(cudaDeviceSynchronize());

        std::jthread t1([&] {
            QubitCode::OneRoundSequential(DMs_Nplus1[1], ORS_DMs_sets[0][1], ORS_DMs_sets[0], DMs_Nplus1[0], d_temps[0], streams_sets[0], handles_sets[0][0], err, errTeleport, QubitCode::correctionMode::noLossW1, 0);
        });

        std::jthread t2([&] {
            QubitCode::OneRoundSequential(DMs_Nplus1[3], ORS_DMs_sets[1][2], ORS_DMs_sets[1], DMs_Nplus1[2], d_temps[1], streams_sets[1], handles_sets[1][0], err, errTeleport, QubitCode::correctionMode::noLossW2, num_anc - 1);

        });

        CHECK_CUDA(cudaDeviceSynchronize());
        // Adding results into (zeroed) input matrix
        combine_2d<<<gridSize_N, blockSize, 0, streams_sets[0][0]>>>(input, ORS_DMs_sets[0][1], ORS_DMs_sets[1][2], DIM);
        // CHECK_CUDA(cudaDeviceSynchronize()); 
    }

    QubitCode::applyAncilla(ORS_DMs_sets[0][0], ancilla_flag[0], DMs_Nplus2[0], DIM, DIM * 4, handles_sets[0][1]);

    CHECK_CUBLAS(cublasSaxpy(handles_sets[0][0], DM_SIZE, &alpha, ORS_DMs_sets[0][0], 1, input, 1));
}