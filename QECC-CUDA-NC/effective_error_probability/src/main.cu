// File to run
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "QubitCode.cuh"
#include "AuxMath.cuh"
#include "kernels.cuh"

#define NUM_STREAMS 4

#define STRCOMPARE(cstring1, cstring2) static_cast<std::string>(cstring1) == cstring2

// So can't have less than 2 rounds in Flag mode
size_t getRounds(int x, const char *mode)
{
    return STRCOMPARE(mode, "flag") ? std::max(static_cast<size_t>(x), static_cast<size_t>(2)) : static_cast<size_t>(x);
}

int main(int argc, char *argv[]) {
    if (argc < 3)
    {
        std::cerr << "Too few arguments! Exiting!\n"
                     "Usage: main.exe [mode (string)] [reencoding error (float)]\n"
                     "\twhere mode can be '1erasure' or 'flag'\n";
        return EXIT_FAILURE;
    }

    try
    {   // MOTE: These paths need to be relative to the executable, not the source code
        readLinesFromFile(corrections, "../../../corrections.txt");                 // Import regular corrections file
        readLinesFromFile(corrections_erasure, "../../../corrections_erasure.txt"); // Import erasure corrections file
        readLinesFromFile(corrections_flag, "../../../corrections_flag.txt");       // Import flag corrections file
        readLinesFromFile(stabilizers, "../../../stabilizers.txt");                 // Import stabilizers file
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    const char *mode = argv[1]; // mode   {"1erasure", "flag"}. NOTE: CL argument at 0 index is name of program itself
    std::istringstream os(argv[2]);
    datatype err;
    os >> err;

    // std::cout << "Here 1\n";

    auto nLinksArray = AuxMath::linspace(1, MAX_DIST, NUMPOINTS);
    const size_t rounds = getRounds(1, mode);               // So 1 round for 1erasure mode, 2 for flag mode. CAN BE CHANGED

    QubitCode::ECInfo *ecInfo = new QubitCode::ECInfo[rounds];
    for (size_t i = 0; i < rounds; i++)
    {
        ecInfo[i].resize(NUMPOINTS);
    }

    // std::cout << "Here 1*\n";


    // PRELIMS DONE
    
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H1[0], DM_H1P));
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H1[1], DM_H1N));
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H2P0_ptr, DM_H2P0));

    // std::cout << "Here 2\n";

    for (int i = 0; i < 2; i++) {
        std::vector<datatype> tmpH1{0.5f, !(i % 2) ? -0.5f : 0.5f, !(i % 2) ? -0.5f : 0.5f, 0.5f};
        // CHECK_CUDA(cudaMemcpyToSymbol(DM_H1[i], tmpH1.data(), tmpH1.size() * sizeof(datatype)));
        CHECK_CUDA(cudaMemcpyToSymbol(i ? DM_H1N : DM_H1P, tmpH1.data(), tmpH1.size() * sizeof(datatype)));
    }

    {
        std::vector<datatype> tmpH2P0{0.5f, 0.5f, 0.0f, 0.0f,
                                      0.5f, 0.5f, 0.0f, 0.0f,
                                      0.0f, 0.0f, 0.0f, 0.0f,
                                      0.0f, 0.0f, 0.0f, 0.0f};
        CHECK_CUDA(cudaMemcpyToSymbol(DM_H2P0, tmpH2P0.data(), tmpH2P0.size() * sizeof(datatype)));
    }

    Matrix I_i(1, 1, {1.0f});
    for (int i = 0; i < NUMQUBITS; i++)
    {
        I_i = kronecker_col_major(I_i, I_h);
    }

    // cuda_global I_d(I_i.mat);       I_d.name = "I_d";
    datatype *I_d;
    CHECK_CUDA(cudaMalloc(&I_d, DM_SIZE * sizeof(datatype)));
    CHECK_CUDA(cudaMemcpy(I_d, I_i.mat.data(), DM_SIZE * sizeof(datatype), cudaMemcpyHostToDevice));

    std::cout << "Here 3\n";


    for (int i = 0; i < 2; i++) {
        Matrix tmp(1, 2, {num, (i % 2) ? -num : num});

        tmp = kronecker_col_major(tmp, I_i);
        // ancilla.emplace_back(tmp.mat);
        CHECK_CUDA(cudaMalloc(&ancilla[i], 2 * DM_SIZE * sizeof(datatype)));
        CHECK_CUDA(cudaMemcpy(ancilla[i], tmp.mat.data(), 2 * DM_SIZE * sizeof(datatype), cudaMemcpyHostToDevice));
    }

    std::cout << "Here 4\n";


    for (int i = 0; i < 4; i++) {
        Matrix tmpH(1, 2, {num, (i % 2) ? -num : num});

        Matrix tmpX(1, 2, {(i < 2) ? 1.0f : 0.0f, (i < 2) ? 0.0f : 1.0f});

        tmpH = kronecker_col_major(tmpH, tmpX);
        tmpH = kronecker_col_major(tmpH, I_i);
        // ancilla_flag.emplace_back(tmpH.mat);
        CHECK_CUDA(cudaMalloc(&ancilla_flag[i], 4 * DM_SIZE * sizeof(datatype)));
        CHECK_CUDA(cudaMemcpy(ancilla_flag[i], tmpH.mat.data(), 4 * DM_SIZE * sizeof(datatype), cudaMemcpyHostToDevice));
    }

    std::cout << "Here 5\n";

    constexpr datatype errTeleFactor = 3;                   // errTeleport = errTeleFactor * eps0.
    constexpr datatype errReencodingFactor = 3;             // eps0 = errReencoding / errReencodingFactor with eps0 being the error of non-teleported 2-qubit gates.
    constexpr double theta = std::numbers::pi / 2;

    auto logicalDM_h = QubitCode::getLogicalDensityMatrix(theta);
    auto logicalKet_h = QubitCode::getLogicalKet(theta);

    // cuda_global logicalDM_d(logicalDM_h);       logicalDM_d.name = "logicalDM_d";
    // cuda_global tempDM_d(logicalDM_h.size());       tempDM_d.name = "tempDM_d";

    datatype *logicalDM_d, *tempDM_d;
    CHECK_CUDA(cudaMalloc(&logicalDM_d, DM_SIZE * sizeof(datatype)));
    CHECK_CUDA(cudaMalloc(&tempDM_d, DM_SIZE * sizeof(datatype)));
    CHECK_CUDA(cudaMemcpy(logicalDM_d, logicalDM_h.data(), DM_SIZE * sizeof(datatype), cudaMemcpyHostToDevice));

    std::cout << "Here 6\n";

    assert(tempDM_d != nullptr && "Failed 1");

    // cuda_global tempX1_d(logicalKet_h.size()), tempY1_d(logicalKet_h.size()), tempZ1_d(logicalKet_h.size());
    // std::vector<std::vector<cuda_global>> ORS_DMs_sets;
    std::vector<std::vector<datatype*>> ORS_DMs_sets(2, std::vector<datatype*>(NUMQUBITS + 2));

    for (auto& DM_set : ORS_DMs_sets) {
        for (auto& DM : DM_set) {
            CHECK_CUDA(cudaMalloc(&DM, DM_SIZE * sizeof(datatype))); 
        }
    }

    std::cout << "Here 7\n";

    assert(tempDM_d != nullptr && "Failed 2");

    // std::cout << "Here 4\n";

    // First one for continuous use, rest are dummies
    std::vector<datatype*> DMs_Nplus1(4);            // NOTE: Very large memory allocations - ensure GPU supports it!

    for (auto& DM : DMs_Nplus1) {
        // DMs_Nplus1.emplace_back(4 * DM_SIZE);
        CHECK_CUDA(cudaMalloc(&DM, 4 * DM_SIZE * sizeof(datatype)));
    }

    // Might want to move it to inside OneRoundSequential (ASYNC alloc - free)
    std::vector<datatype*> DMs_Nplus2(2);            // NOTE: Very large memory allocations - ensure GPU supports it!

    for (auto& DM : DMs_Nplus2) {
        // DMs_Nplus2.emplace_back(16 * DM_SIZE);
        CHECK_CUDA(cudaMalloc(&DM, 16 * DM_SIZE * sizeof(datatype)));
    }

    // datatype** d_temps_1, **d_temps_2;         // Had to do this since kernels can't take std::vector<cuda_global> as argument
    // std::pair<datatype**, datatype**> d_temps = std::make_pair(d_temps_1, d_temps_2);
    std::array<datatype**, 2> d_temps;

    std::cout << "Here 8\n";

    for (auto d_temp: d_temps) {                        // This looks VERY sus... we're not mallocing to d_temp[i]
        cudaMalloc(&d_temp, 12 * sizeof(datatype*));
        for (int i = 0; i < 12; ++i) {
            datatype* ptr;
            CHECK_CUDA(cudaMalloc(&ptr, (DIM * 4) * (DIM * 4) * sizeof(datatype)));
            CHECK_CUDA(cudaMemcpy(&d_temp[i], &ptr, sizeof(datatype*), cudaMemcpyHostToDevice));
        }
    }

    std::array<std::vector<cudaStream_t>, 2> streams_sets({
        std::vector<cudaStream_t>(NUM_STREAMS),
        std::vector<cudaStream_t>(NUM_STREAMS)
    });

    for (auto &streams : streams_sets) {
        for (auto &stream : streams) {
            CHECK_CUDA(cudaStreamCreate(&stream));
        }
    }

    std::array<std::vector<cublasHandle_t>, 2> handles_sets({
        std::vector<cublasHandle_t>(NUM_STREAMS),
        std::vector<cublasHandle_t>(NUM_STREAMS)
    });


    for (auto &handles : handles_sets) {
        for (auto &handle : handles) {
            CHECK_CUBLAS(cublasCreate(&handle));
        }
    }

    for (unsigned i = 0; i < 2; i++) {
        for (unsigned j = 0; j < NUM_STREAMS; j++) {
            CHECK_CUBLAS(cublasSetStream(handles_sets[i][j], streams_sets[i][j]));
        }
    }
    // May combine these into a single loop later

    // std::cout << "Here 6\n";

    datatype etrans;
    std::vector<datatype*> traces_d(rounds), fidels_d(rounds);

    for (int i = 0; i < rounds; i++) {
        CHECK_CUDA(cudaMalloc(&traces_d[i], NUMPOINTS * sizeof(datatype)));
        CHECK_CUDA(cudaMalloc(&fidels_d[i], NUMPOINTS * sizeof(datatype)));
    }

    dim3 gridSize = getGridSize(DIM, DIM);
    
    CHECK_CUDA(cudaDeviceSynchronize());        // Just for safety

    assert(tempDM_d != nullptr && "Failed 3");


    for (unsigned i = 0; i < nLinksArray.size(); i++) {
        // logicalDM_d.size--;
        assert(tempDM_d != nullptr && "Failed 4");   
        // tempDM_d = logicalDM_d;
        CHECK_CUDA(cudaMemcpy(tempDM_d, logicalDM_d, DM_SIZE * sizeof(datatype), cudaMemcpyDeviceToDevice));
        // assert(tempDM_d.size == DM_SIZE && "Not matched");
        assert(tempDM_d != nullptr && "Failed 5");   


        for (unsigned j = 0; j < rounds; j++) {
            etrans = 1 - std::pow(1 - err, nLinksArray[i]) * (1 - err / errReencodingFactor);

            for (unsigned k = 0; k < NUMQUBITS; k++) {
                // Single-Qubit Uniform Depolarization
                std::cout << "k = " << k << std::endl;
                assert(tempDM_d != nullptr && "Failed 6");   
                // CHECK_CUDA(cudaGetLastError());

                // These run concurrently
                X_dm<<<gridSize, blockSize, 0, streams_sets[0][1]>>>(d_temps[0][0], tempDM_d, k, DIM);
                // CHECK_CUDA(cudaStreamSynchronize(streams_sets[0][1]));
                // CHECK_CUDA(cudaGetLastError());
                Y_dm<<<gridSize, blockSize, 0, streams_sets[0][2]>>>(d_temps[0][1], tempDM_d, k, DIM);
                // CHECK_CUDA(cudaStreamSynchronize(streams_sets[0][2]));
                // CHECK_CUDA(cudaGetLastError());
                Z_dm<<<gridSize, blockSize, 0, streams_sets[0][3]>>>(d_temps[0][2], tempDM_d, k, DIM);
                // CHECK_CUDA(cudaStreamSynchronize(streams_sets[0][3]));
                // CHECK_CUDA(cudaGetLastError());
                CHECK_CUDA(cudaDeviceSynchronize());

                combine_and_scale<<<gridSize, blockSize, 0, streams_sets[0][0]>>>(tempDM_d, d_temps[0], 3, 1 - etrans, etrans / 3, DIM);
                CHECK_CUDA(cudaDeviceSynchronize());
            }

            // CHECK_CUDA(cudaMalloc(&d_temps, d_temps_size));
            // This was done to save memory at all other places, but seems unnecessary unless VRAM is overloaded
            if (STRCOMPARE(mode, "flag")) {
                QubitCode::OneRoundSequential_flag(tempDM_d, ORS_DMs_sets, DMs_Nplus2, DMs_Nplus1, d_temps, streams_sets, handles_sets, 
                                                   err / errReencodingFactor, 
                                                   errTeleFactor * err / errReencodingFactor, 
                                                   errTeleFactor * err / errReencodingFactor);
            }
            else {
                QubitCode::OneRoundSequential(DMs_Nplus1[1], tempDM_d, ORS_DMs_sets[0], DMs_Nplus1[0], d_temps[0], streams_sets[0], handles_sets[0][0], 
                                              err / errReencodingFactor, 
                                              errTeleFactor * err / errReencodingFactor, 
                                              QubitCode::correctionMode::oneLoss, 0);
            }
            // CHECK_CUDA(cudaFree(d_temps));

            CHECK_CUDA(cudaDeviceSynchronize());

            traceDM(handles_sets[0][1], DIM, tempDM_d, I_d, traces_d[j], i);
            fidelDM(handles_sets[0][2], DIM, tempDM_d, logicalDM_d, fidels_d[j], i);

            // trace_dm<<<gridSize, blockSize, 0, streams[1]>>>(traces_d[j].get(), correctedDM_d.get(), DIM);
            // fidel_dm<<<gridSize, blockSize, 0, streams[2]>>>(fidels_d[j].get(), correctedDM_d.get(), logicalDM_d.get(), DIM);

            // tempDM_d = correctedDM_d;       // This could be made async too
            // CHECK_CUDA(cudaMemcpyAsync(tempDM_d.get(), correctedDM_d.get(), logicalDM_h.size() * sizeof(datatype), cudaMemcpyDeviceToDevice, streams_sets[0][3]));
            
            CHECK_CUDA(cudaDeviceSynchronize());

        }
    }

    for (auto d_temp: d_temps) {
        datatype* h_ptrs[12];
        CHECK_CUDA(cudaMemcpy(h_ptrs, d_temp, 12 * sizeof(datatype*), cudaMemcpyDeviceToHost));     // ERROR

        for (auto* ptr: h_ptrs) {
            CHECK_CUDA(cudaFree(ptr));  // Free each vector
        }

        CHECK_CUDA(cudaFree(d_temp));
    }


    for (unsigned j = 0; j < rounds; j++) {
        CHECK_CUDA(cudaMemcpyAsync(ecInfo[j].grandTraces.data(), traces_d[j], NUMPOINTS * sizeof(datatype), cudaMemcpyDeviceToHost, streams_sets[0][1]));
        CHECK_CUDA(cudaMemcpyAsync(ecInfo[j].grandErrors.data(), fidels_d[j], NUMPOINTS * sizeof(datatype), cudaMemcpyDeviceToHost, streams_sets[0][2]));
        std::transform(ecInfo[j].grandErrors.begin(), ecInfo[j].grandErrors.end(), ecInfo[j].grandErrors.begin(), [](datatype x) { return 1 - x; });
    }
    
    CHECK_CUDA(cudaDeviceSynchronize());

    // FREE EVERYTHING !!!
    CHECK_CUDA(cudaFree(I_d));
    CHECK_CUDA(cudaFree(logicalDM_d));
    CHECK_CUDA(cudaFree(tempDM_d));

    for (auto& DM_set : ORS_DMs_sets) {
        for (auto& DM : DM_set) {
            CHECK_CUDA(cudaFree(DM));
        }
    }

    for (auto &handles : handles_sets) {
        for (auto &handle : handles) {
            CHECK_CUBLAS(cublasDestroy(handle));
        }
    }

    for (auto &streams : streams_sets) {
        for (auto &stream : streams) {
            CHECK_CUDA(cudaStreamDestroy(stream));
        }
    }

    for (auto& DM : DMs_Nplus1) {
        // DMs_Nplus1.emplace_back(4 * DM_SIZE);
        CHECK_CUDA(cudaFree(DM));
    }

    for (auto& DM : DMs_Nplus2) {
        // DMs_Nplus1.emplace_back(4 * DM_SIZE);
        CHECK_CUDA(cudaFree(DM));
    }

    for (int i = 0; i < rounds; i++) {
        CHECK_CUDA(cudaFree(traces_d[i]));
        CHECK_CUDA(cudaFree(fidels_d[i]));
    }

    delete[] ecInfo;

    return EXIT_SUCCESS;
}