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

size_t getRounds(int x, const char *mode);
std::string datetime();
std::string to_string_with_precision(const double a_value, const int n);

int main(int argc, char *argv[]) {
    if (argc < 3)
    {
        std::cerr << "Too few arguments! Exiting!\n"
                     "Usage: main.exe [mode (string)] [reencoding error (float)]\n"
                     "\twhere mode can be '1erasure' or 'flag'\n";
        return EXIT_FAILURE;
    }

    std::string project_root = "D:/QEC_Codes/QECC-CUDA/";

    try
    {   // MOTE: These paths need to be relative to the executable, not the source code
        readLinesFromFile(corrections, project_root + "corrections.txt");                 // Import regular corrections file
        readLinesFromFile(corrections_erasure, project_root + "corrections_erasure.txt"); // Import erasure corrections file
        readLinesFromFile(corrections_flag, project_root + "corrections_flag.txt");       // Import flag corrections file
        readLinesFromFile(stabilizers, project_root + "stabilizers.txt");                 // Import stabilizers file
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

    const auto t1 = std::chrono::high_resolution_clock::now();
    auto nLinksArray = AuxMath::linspace(1, MAX_DIST, NUMPOINTS);
    const size_t rounds = getRounds(1, mode);               // So 1 round for 1erasure mode, 2 for flag mode. CAN BE CHANGED

    QubitCode::ECInfo *ecInfo = new QubitCode::ECInfo[rounds];
    for (size_t i = 0; i < rounds; i++)
    {
        ecInfo[i].resize(NUMPOINTS);
    }

    // PRELIMS DONE
    
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H1[0], DM_H1P));
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H1[1], DM_H1N));
    // CHECK_CUDA(cudaGetSymbolAddress((void**)DM_H2P0_ptr, DM_H2P0));

    // for (int i = 0; i < 2; i++) {
    //     std::vector<datatype> tmpH1{0.5f, !(i % 2) ? -0.5f : 0.5f, !(i % 2) ? -0.5f : 0.5f, 0.5f};
    //     // CHECK_CUDA(cudaMemcpyToSymbol(DM_H1[i], tmpH1.data(), tmpH1.size() * sizeof(datatype)));
    //     CHECK_CUDA(cudaMemcpyToSymbol(i ? DM_H1N : DM_H1P, tmpH1.data(), tmpH1.size() * sizeof(datatype)));
    // }

    std::vector<datatype> tmpH1P{0.5f, 0.5f, 0.5f, 0.5f};
    DM_H1P = std::move(cuda_global(tmpH1P));

    
    std::vector<datatype> tmpH2P0{0.5f, 0.5f, 0.0f, 0.0f,
                                  0.5f, 0.5f, 0.0f, 0.0f,
                                  0.0f, 0.0f, 0.0f, 0.0f,
                                  0.0f, 0.0f, 0.0f, 0.0f};
    // CHECK_CUDA(cudaMemcpyToSymbol(DM_H2P0, tmpH2P0.data(), tmpH2P0.size() * sizeof(datatype)));
    DM_H2P0 = std::move(cuda_global(tmpH2P0));
    

    Matrix I_i(1, 1, {1.0f});
    for (int i = 0; i < NUMQUBITS; i++)
    {
        I_i = kronecker_col_major(I_i, I_h);
    }
    cuda_global I_d(I_i.mat);


    for (int i = 0; i < 2; i++) {
        Matrix tmp(1, 2, {num, (i % 2) ? -num : num});

        tmp = kronecker_col_major(tmp, I_i);
        ancilla.emplace_back(tmp.mat);
    }


    for (int i = 0; i < 4; i++) {
        Matrix tmpH(1, 2, {num, (i % 2) ? -num : num});

        Matrix tmpX(1, 2, {(i < 2) ? 1.0f : 0.0f, (i < 2) ? 0.0f : 1.0f});

        tmpH = kronecker_col_major(tmpH, tmpX);
        tmpH = kronecker_col_major(tmpH, I_i);
        ancilla_flag.emplace_back(tmpH.mat);
    }

    constexpr datatype errTeleFactor = 3;                   // errTeleport = errTeleFactor * eps0.
    constexpr datatype errReencodingFactor = 3;             // eps0 = errReencoding / errReencodingFactor with eps0 being the error of non-teleported 2-qubit gates.
    constexpr double theta = std::numbers::pi / 2;

    const std::string dataPath = project_root + "binary_data/"; // Output Binary Data created here
    std::filesystem::create_directory(dataPath);
    const std::string dataSpecificPath = dataPath + mode + (std::string) "_err" + to_string_with_precision(err, 0) + (std::string) "_reencodingFactor" + std::to_string((int)errReencodingFactor) + (std::string) "_teleFactor" + std::to_string((int)errTeleFactor) + (std::string) "/";
    std::filesystem::create_directory(dataSpecificPath);


    auto logicalDM_h = QubitCode::getLogicalDensityMatrix(theta);
    auto logicalKet_h = QubitCode::getLogicalKet(theta);

    cuda_global logicalDM_d(logicalDM_h);
    cuda_global tempDM_d(logicalDM_h.size());

    // cuda_global tempX1_d(logicalKet_h.size()), tempY1_d(logicalKet_h.size()), tempZ1_d(logicalKet_h.size());
    std::vector<std::vector<cuda_global>> ORS_DMs_sets;

    for (int i = 0; i < 2; i++) {
        std::vector<cuda_global> vec_d;
        for (int j = 0; j < NUMQUBITS + 2; j++) {
            vec_d.emplace_back(logicalDM_h.size());
        }
        ORS_DMs_sets.push_back(std::move(vec_d));
    }

    // First one for continuous use, rest are dummies
    std::vector<cuda_global> DMs_Nplus1;            // NOTE: Very large memory allocations - ensure GPU supports it!

    for (int j = 0; j < 4; j++) {
        DMs_Nplus1.emplace_back(4 * DM_SIZE);
    }

    // Might want to move it to inside OneRoundSequential (ASYNC alloc - free)
    std::vector<cuda_global> DMs_Nplus2;

    for (int j = 0; j < 2; j++) {
        DMs_Nplus2.emplace_back(16 * DM_SIZE);
    }

    // datatype** d_temps_1, **d_temps_2;         // Had to do this since kernels can't take std::vector<cuda_global> as argument
    // std::pair<datatype**, datatype**> d_temps = std::make_pair(d_temps_1, d_temps_2);

    std::vector<std::vector<cuda_global>> d_temps(2);

    for (auto& d_temp: d_temps) {
        for (int i = 0; i < 12; ++i) {
            // datatype* ptr;
            // cudaMalloc(&ptr, (DIM * 4) * (DIM * 4) * sizeof(datatype));
            // cudaMemcpy(&d_temp[i], &ptr, sizeof(datatype*), cudaMemcpyHostToDevice);
            // CHECK_CUDA(cudaMalloc(&d_temp[i], (DIM * 4) * (DIM * 4) * sizeof(datatype)));
            d_temp.emplace_back(DM_SIZE);
        }
    }

    std::vector<std::vector<cudaStream_t>> streams_sets(2, std::vector<cudaStream_t>(NUM_STREAMS));

    for (auto &streams : streams_sets) {
        for (auto &stream : streams) {
            CHECK_CUDA(cudaStreamCreate(&stream));
        }
    }

    std::vector<std::vector<cublasHandle_t>> handles_sets(2, std::vector<cublasHandle_t>(NUM_STREAMS));


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

    datatype etrans;
    std::vector<cuda_global> traces_d, fidels_d;

    for (int i = 0; i < rounds; i++) {
        traces_d.emplace_back(NUMPOINTS);
        fidels_d.emplace_back(NUMPOINTS);
    }

    dim3 gridSize = getGridSize(DIM, DIM);
    
    CHECK_CUDA(cudaDeviceSynchronize());        // Just for safety


    for (unsigned i = 0; i < nLinksArray.size(); i++) {
        //std::cout << "i = " << i << std::endl;

        tempDM_d = logicalDM_d;

        //assert(tempDM_d.size == DM_SIZE && "Not matched");
        //assert(tempDM_d.get() != nullptr && "Failed");

        CHECK_CUDA(cudaDeviceSynchronize());

        for (unsigned j = 0; j < rounds; j++) {
            etrans = 1 - std::pow(1 - err, nLinksArray[i]) * (1 - err / errReencodingFactor);

            for (unsigned k = 0; k < NUMQUBITS; k++) {
                // Single-Qubit Uniform Depolarization
                // std::cout << "k = " << k << std::endl;
                assert(tempDM_d.get() != nullptr && "Failed 6");   
                // CHECK_CUDA(cudaGetLastError());

                // These run concurrently
                X_dm<<<gridSize, blockSize, 0, streams_sets[0][1]>>>(d_temps[0][0].get(), tempDM_d.get(), k, DIM);

                Y_dm<<<gridSize, blockSize, 0, streams_sets[0][2]>>>(d_temps[0][1].get(), tempDM_d.get(), k, DIM);

                Z_dm<<<gridSize, blockSize, 0, streams_sets[0][3]>>>(d_temps[0][2].get(), tempDM_d.get(), k, DIM);

                CHECK_CUDA(cudaDeviceSynchronize());
    

                combine_and_scale_3<<<gridSize, blockSize, 0, streams_sets[0][0]>>>(tempDM_d.get(), d_temps[0][0].get(), d_temps[0][1].get(), d_temps[0][2].get(), 1 - etrans, etrans / 3, DIM);
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
            
            CHECK_CUDA(cudaDeviceSynchronize());

        }
    }

    // for (auto d_temp: d_temps) {
    //     datatype* h_ptrs[12];
    //     CHECK_CUDA(cudaMemcpy(h_ptrs, d_temp, 12 * sizeof(datatype*), cudaMemcpyDeviceToHost));

    //     for (auto* ptr: h_ptrs) {
    //         CHECK_CUDA(cudaFree(ptr));  // Free each vector
    //     }

    //     CHECK_CUDA(cudaFree(d_temp));
    // }

    // for (auto d_temp: d_temps) {
    //     for (auto i = 0; i < 12; i++) {
    //         CHECK_CUDA(cudaFree(d_temp[i]));  // Free each vector
    //     }

    //     CHECK_CUDA(cudaFree(d_temp));
    // }


    for (unsigned j = 0; j < rounds; j++) {
        CHECK_CUDA(cudaMemcpyAsync(ecInfo[j].grandTraces.data(), traces_d[j].get(), NUMPOINTS * sizeof(datatype), cudaMemcpyDeviceToHost, streams_sets[0][1]));
        CHECK_CUDA(cudaMemcpyAsync(ecInfo[j].grandErrors.data(), fidels_d[j].get(), NUMPOINTS * sizeof(datatype), cudaMemcpyDeviceToHost, streams_sets[0][2]));
        std::transform(ecInfo[j].grandErrors.begin(), ecInfo[j].grandErrors.end(), ecInfo[j].grandErrors.begin(), [](const datatype& x) { return 1 - x; });
    }
    
    CHECK_CUDA(cudaDeviceSynchronize());

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

    const auto instantaneous_datetime = datetime();
    for (size_t i = 0; i < rounds; i++)
    {
        // std::cout << "ecInfo loop" << std::endl;
        ecInfo[i].print();
        ecInfo[i].save(dataSpecificPath + (std::string) "ecInfo" + std::to_string(i) + (std::string) "_" + instantaneous_datetime + (std::string) ".dat");
    }

    // This binary file can be read in Wolfram Mathematica by writing BinaryReadList["err.dat", "Real64"] and then you can proceed to plot it.
    std::ofstream outfile_nLinksArray((dataSpecificPath + (std::string) "nLinksArray_" + instantaneous_datetime + (std::string) ".dat"), std::ios::binary);
    outfile_nLinksArray.write((char *)nLinksArray.data(), NUMPOINTS * sizeof(double));
    outfile_nLinksArray.close();

    delete[] ecInfo;

    return EXIT_SUCCESS;
}

std::string datetime()
{
    time_t rawtime;
    struct tm *timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y-%m-%d %H-%M-%S", timeinfo);
    return std::string(buffer);
}

// Converts float/double to string with the desired precision
std::string to_string_with_precision(const double a_value, const int n)
{
    // Reference: https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
    std::ostringstream out;
    out.precision(n);
    out << std::scientific << a_value;
    return out.str();
}

// So can't have less than 2 rounds in Flag mode
size_t getRounds(int x, const char *mode)
{
    return STRCOMPARE(mode, "flag") ? std::max(static_cast<size_t>(x), static_cast<size_t>(2)) : static_cast<size_t>(x);
}