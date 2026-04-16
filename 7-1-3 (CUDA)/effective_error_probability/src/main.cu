#include <array>
#include <chrono> // https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream> // Provides std::ofstream
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>
// #include <iterator>
// #include <algorithm>

#include "QubitCode.cuh"
#include "AuxMath.cuh"

__constant__ datatype I_d[4], X_d[4], Y_d[4], Z_d[4];       // COLUMN MAJOR!!!

#define STRCOMPARE(cstring1, cstring2) static_cast<std::string>(cstring1) == cstring2


size_t getRounds(int x, const char *mode);
std::string datetime();
std::string to_string_with_precision(const double a_value, const int n);

void readLinesFromFile(std::vector<std::string> &lines, const std::string &filename)    // Modified to not make a large copy of lines when returning
{                                                                                       // Please modifiy similar cases/functions
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::copy(std::istream_iterator<std::string>(file), std::istream_iterator<std::string>(), std::back_inserter(lines));
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Too few arguments! Exiting!\n"
                     "Usage: main.exe [mode (string)] [reencoding error (float)]\n"
                     "\twhere mode can be '1erasure' or 'flag'\n";
        return EXIT_FAILURE;
    }

    try
    {
        readLinesFromFile(corrections, "corrections.txt");                 // Import regular corrections file
        readLinesFromFile(corrections_erasure, "corrections_erasure.txt"); // Import erasure corrections file
        readLinesFromFile(corrections_flag, "corrections_flag.txt");       // Import flag corrections file
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    const char *mode = argv[1]; // mode ∈ {"1erasure", "flag"}. NOTE: CL argument at 0 index is name of program itself
    std::istringstream os(argv[2]);
    double err;
    os >> err;

    if (!(STRCOMPARE(mode, "1erasure") || STRCOMPARE(mode, "flag")))
    {
        std::cerr << "Usage: main.exe [mode (string)] [reencoding error (float)]\n";
        std::cerr << "\twhere mode can be '1erasure' or 'flag'\n";
        return EXIT_FAILURE;
    }
    if (argc > 3)
    {
        std::cerr << "Warning: Extra arguments ignored\n";
    }

    QubitCode::Form_Satbilizers("stabilizers.txt");
    auto cArrLen = FastPauli::CArrLen<NUMQUBITS>;

    // NOTE: Eigen's .data() passes data as ColMajor
    CHECK_CUDA(cudaMemcpyToSymbol(I_d, I_h.data(), 4 * sizeof(datatype)));
    CHECK_CUDA(cudaMemcpyToSymbol(X_d, X_h.data(), 4 * sizeof(datatype)));
    CHECK_CUDA(cudaMemcpyToSymbol(Y_d, Y_h.data(), 4 * sizeof(datatype)));
    CHECK_CUDA(cudaMemcpyToSymbol(Z_d, Z_h.data(), 4 * sizeof(datatype)));


    std::vector<datatype *> XYZ_d(NUMQUBITS);

    Operator_t tmp_XYZ_h;
    for (int i = 0; i < NUMQUBITS; i++)
    {
        Eigen::Matrix<datatype, dynamic, dynamic> X_i = {1}, Y_i = {1}, Z_i = {1};
        for (int j = 0; j < NUMQUBITS; j++) {
            X_i = kroneckerProduct(X_i, (i == j) ? X_h : I_h).eval();
            Y_i = kroneckerProduct(Y_i, (i == j) ? Y_h : I_h).eval();
            Z_i = kroneckerProduct(Z_i, (i == j) ? Z_h : I_h).eval();
        }

        tmp_XYZ_h = X_i + Y_i + Z_i;
        CHECK_CUDA(cudaMalloc((void **)&XYZ_d[i], cArrLen * sizeof(datatype)));
        copyUpperTriangularColumnMajor(tmp_XYZ_h, XYZ_d[i], cArrLen);
    }


    const auto t1 = std::chrono::high_resolution_clock::now();
    constexpr double errTeleFactor = 3;       // errTeleport = errTeleFactor * eps0.
    constexpr double errReencodingFactor = 3; // eps0 = errReencoding / errReencodingFactor with eps0 being the error of non-teleported 2-qubit gates.
    constexpr double theta = std::numbers::pi / 2;
    auto nLinksArray = AuxMath::linspace_host(1, MAX_DIST, NUMPOINTS);
    const size_t rounds = getRounds(1, mode); // So 1 round for 1erasure mode, 2 for flag mode. CAN BE CHANGED

    std::cout << "  \x1b[0;94mmode: " << mode << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94mNUMPOINTS: " << NUMPOINTS << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94merr: " << to_string_with_precision(err, 0) << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94merrTeleFactor: " << errTeleFactor << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94merrReencodingFactor: " << errReencodingFactor << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94mtheta: " << theta << "\x1b[0m" << std::endl;
    std::cout << "  \x1b[0;94mrounds: " << rounds << "\x1b[0m" << std::endl;

    const std::string dataPath = "./binary_data/"; // Output Binary Data created here
    std::filesystem::create_directory(dataPath);
    const std::string dataSpecificPath = dataPath + mode + (std::string) "_err" + to_string_with_precision(err, 0) + (std::string) "_reencodingFactor" + std::to_string((int)errReencodingFactor) + (std::string) "_teleFactor" + std::to_string((int)errTeleFactor) + (std::string) "/";
    std::filesystem::create_directory(dataSpecificPath);

    QubitCode::ECInfo *ecInfo = new QubitCode::ECInfo[rounds];
    for (size_t i = 0; i < rounds; i++)
    {
        ecInfo[i].resize(NUMPOINTS);
    }

    // CHECK DONE

    // Logical |+> is used as inital input. Fidelity is also calculated wrt logicalcArr below
    std::vector<datatype> logicalcArr = QubitCode::getLogicalDensityMatrix(theta); // cArr of |+>
    std::vector<datatype> logicalKet = QubitCode::getLogicalKet(theta);
    std::vector<datatype> tempcArrCodeword(FastPauli::CArrLen<NUMQUBITS>, 0.0);

    std::cout << "fidelity = " << FastPauli::fidel_r(logicalKet, logicalcArr, NUMQUBITS) << std::endl;
    std::cout << "trace = " << FastPauli::trace_r(logicalcArr, NUMQUBITS) << std::endl;

    std::vector<datatype> correctedcArr(FastPauli::CArrLen<NUMQUBITS>, 0.0);
    datatype etrans;

    int i, j;

    for (i = 0; i < nLinksArray.size(); i++)
    {
        // std::cout << "Start of i for loop" << std::endl;
        tempcArrCodeword = logicalcArr;
        for (j = 0; j < rounds; j++)
        {
            etrans = 1 - std::pow(1 - err, nLinksArray[i]) * (1 - err / errReencodingFactor);

            // "Uniform" depolarisation on each data qubit applied
            for (int k = 0; k < NUMQUBITS; k++)
            {
                NoisyPauli::SingleQubitDepolarise_r(tempcArrCodeword, NUMQUBITS, k, etrans);
            }
            // std::cout << "Inside MAIN flag loop (first)" << std::endl;
            correctedcArr = QubitCode::OneRoundSequential_flag(logicalcArr,
                                                                    tempcArrCodeword,
                                                                    err / errReencodingFactor,
                                                                    errTeleFactor * err / errReencodingFactor,
                                                                    errTeleFactor * err / errReencodingFactor);
            // std::cout << "Inside MAIN flag loop (last)" << std::endl;
            ecInfo[j].grandErrors[i] = (1 - FastPauli::fidel_r(logicalKet, correctedcArr, NUMQUBITS));
            ecInfo[j].grandTraces[i] = FastPauli::trace_r(correctedcArr, NUMQUBITS);

            tempcArrCodeword = correctedcArr;
        }
    }

    const auto instantaneous_datetime = datetime();
    for (size_t i = 0; i < rounds; i++)
    {
        // std::cout << "ecInfo loop" << std::endl;
        ecInfo[i].save(dataSpecificPath + (std::string) "ecInfo" + std::to_string(i) + (std::string) "_" + instantaneous_datetime + (std::string) ".dat");
    }

    // This binary file can be read in Wolfram Mathematica by writing BinaryReadList["err.dat", "Real64"] and then you can proceed to plot it.
    std::ofstream outfile_nLinksArray((dataSpecificPath + (std::string) "nLinksArray_" + instantaneous_datetime + (std::string) ".dat"), std::ios::binary);
    outfile_nLinksArray.write((char *)nLinksArray.data(), NUMPOINTS * sizeof(double));
    outfile_nLinksArray.close();

    // std::cout << "After nLinksArray" << std::endl;

    delete[] ecInfo;
    for (auto x : XYZ_d)
    {
        CHECK_CUDA(cudaFree(x));
    }

    const auto t2 = std::chrono::high_resolution_clock::now();
    const auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "\n\033[1;34mElapsed time: \033[0m\033[1;4;33m" << ms_int.count() << "ms\033[0m\n";
    std::cout << "instantaneous_datetime: " << instantaneous_datetime << std::endl;
    std::cout << '\a' << std::endl;
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