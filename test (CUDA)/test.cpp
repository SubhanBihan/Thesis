// Just a file to test basic code

#include <iostream>
#include <filesystem>
#include <string>
#include <format>
#include <vector>
#include <cassert>
#include <random>

std::string findFileWithPartName(const std::string &directoryPath, const std::string &partOfFileName)
{
    namespace fs = std::filesystem;

    // Iterate through the directory
    for (const auto &entry : fs::directory_iterator(directoryPath))
    {
        // Check if the entry is a file
        if (entry.is_regular_file())
        {
            // Get the filename as a string
            std::string filename = entry.path().filename().string();

            // Check if the filename contains the partOfFileName
            if (filename.find(partOfFileName) != std::string::npos)
            {
                // Return the full path of the file
                return entry.path().string();
            }
        }
    }

    // Return an empty string if no file is found
    return "";
}

void apply_permutation(std::vector<double> &v, int *indices, const int vSize)
{ // Reference: https://stackoverflow.com/questions/67751784/how-to-do-in-place-sorting-a-list-according-to-a-given-index-in-c.
    for (int i = 0; i < vSize; i++)
    {
        auto current = i;
        while (i != indices[current])
        {
            auto next = indices[current];
            std::swap(v[current], v[next]);
            indices[current] = current;
            current = next;
        }
        indices[current] = current;
    }
}

int fastMod(const int i, const int t)
{
    return (i) & ((1 << t) - 1);
}

void ketKron(std::vector<double> &result, const std::vector<double> &ket1, const std::vector<double> &ket2, const int numQubits1, const int numQubits2)
{
    const int numQubitsTot = numQubits1 + numQubits2;
    const int ntot = 1 << numQubitsTot;
    // const int n1 = 1 << numQubits1; // Not used.
    const int n2 = 1 << numQubits2;
    for (int i = 0; i < ntot; i++)
    {
        result[i] = ket1[i / n2] * ket2[fastMod(i, numQubits2)];
    }
}

void ketX_r(std::vector<double> &ket, const int numQubits, const int targetQubit)
{
    int i;
    const int n = 1 << numQubits; // Matrix dimension.
    const int tmp0 = (n >> (targetQubit + 1));
    int *ketIndex = new int[n];
    for (i = 0; i < n; i++)
    {
        ketIndex[i] = i ^ tmp0;
    }
    apply_permutation(ket, ketIndex, n);
    delete[] ketIndex;
}

void ketZ_r(std::vector<double>& ket, const int numQubits, const int targetQubit) {
    int i, ketIndex;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int n2 = n >> 1;
    const int tmp0 = (n >> (targetQubit + 1));
    for (i = 0; i < n2; i++) {
        ketIndex = n + 1 - (2 * (i + 1) - tmp0 + fastMod(-i - 1, numQubits - targetQubit - 1)) - 1;
        ket[ketIndex] = -ket[ketIndex];
    }
}

void ketY_r(std::vector<double>& ket, const int numQubits, const int targetQubit) {  // Up to some global phase!
    ketZ_r(ket, numQubits, targetQubit);
    ketX_r(ket, numQubits, targetQubit);
}

void ketH_r(std::vector<double> &ket, const int numQubits, const int targetQubit)
{
    const int n = 1 << numQubits; // Total size of state vector (2^numQubits)
    // For leftmost indexing, shift by (numQubits - 1 - targetQubit)
    const int targetPower = 1 << (numQubits - 1 - targetQubit);
    const double normFactor = 1.0 / std::numbers::sqrt2; // 1/√2 normalization factor

    // Process pairs of elements that differ only in targetQubit
    for (int i = 0; i < n; i++)
    {
        // Skip if targetQubit bit is set to ensure we process each pair only once
        if (i & targetPower)
            continue;

        // j is i with targetQubit bit flipped
        const int j = i | targetPower;

        // Store original values
        const double vi = ket[i];
        const double vj = ket[j];

        // Apply Hadamard transformation:
        // |0⟩ -> (|0⟩ + |1⟩)/√2
        // |1⟩ -> (|0⟩ - |1⟩)/√2
        ket[i] = (vi + vj) * normFactor;
        ket[j] = (vi - vj) * normFactor;
    }
}

static void ketCZ_r(std::vector<double> &ket, const int numQubits, const int controlQubit, const int targetQubit)
{
    const int n = 1 << numQubits; // Total size of state vector (2^numQubits)

    // Calculate bit positions (using leftmost indexing convention)
    const int controlPower = 1 << (numQubits - 1 - controlQubit);
    const int targetPower = 1 << (numQubits - 1 - targetQubit);

    // Iterate through all basis states
    for (int i = 0; i < n; i++)
    {
        // Check if both control and target qubits are 1
        if ((i & controlPower) && (i & targetPower))
        {
            // Apply phase flip (-1) when both qubits are 1
            ket[i] = -ket[i];
        }
    }
}

void ketCX_r(std::vector<double> &ket, const int numQubits, const int controlQubit, const int targetQubit)
{
    const int n = 1 << numQubits; // Total size of state vector (2^numQubits)

    // Calculate bit positions (using leftmost indexing convention)
    const int controlPower = 1 << (numQubits - 1 - controlQubit);
    const int targetPower = 1 << (numQubits - 1 - targetQubit);

    // Iterate through all basis states
    for (int i = 0; i < n; i++)
    {
        // Check if control qubit is 1 and we haven't processed this pair yet
        if ((i & controlPower) && !(i & targetPower))
        {
            // j is i with targetQubit flipped
            const int j = i | targetPower;

            // Swap amplitudes when control is 1
            std::swap(ket[i], ket[j]);
        }
    }
}

static void applyAncillaKet(std::vector<double> &out_ket, const int numAncillaQubits, const int numQubits, const std::vector<double> &ancillaKet, const std::vector<double> &in_ket)
{
    // Validate input sizes
    const int ancillaSize = 1 << numAncillaQubits;
    const int mainSize = 1 << numQubits;
    const int outSize = 1 << (numQubits - numAncillaQubits);

    assert(ancillaKet.size() == ancillaSize);
    assert(in_ket.size() == mainSize);

    // Resize output vector
    out_ket.resize(outSize);

    // Zero out the output vector
    std::fill(out_ket.begin(), out_ket.end(), 0.0);

    // Iterate through all possible main (non-ancilla) states
    for (int reducedMainIdx = 0; reducedMainIdx < outSize; ++reducedMainIdx)
    {
        // Find full state index by shifting reduced index into correct bit position
        for (int ancillaIdx = 0; ancillaIdx < ancillaSize; ++ancillaIdx)
        {
            // Check if ancilla amplitude is significant
            if (std::abs(ancillaKet[ancillaIdx]) < 1e-10)
                continue;

            // Construct full state index by inserting ancilla bits
            const int fullIdx = (ancillaIdx << (numQubits - numAncillaQubits)) | reducedMainIdx;

            // Accumulate amplitude
            out_ket[reducedMainIdx] += ancillaKet[ancillaIdx] * in_ket[fullIdx];
        }
    }
}

int simulateMultipleEvents(const std::vector<double> &eventProbabilities)
{
    // Create a random device and engine
    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate a random number between 0 and 1
    double randomValue = dist(engine);

    // Iterate through cumulative probabilities to find the event
    double cumulativeProbability = 0.0;
    for (size_t i = 0; i < eventProbabilities.size(); ++i)
    {
        cumulativeProbability += eventProbabilities[i];
        if (randomValue <= cumulativeProbability)
        {
            return i + 1; // Return event code (1-based index)
        }
    }
    return 0; // No event occurred
}

double trace_ket(const std::vector<double> &ket, const int numQubits)
{
    // std::vector<double> cArr(numQubits);
    // FastPauli::ketToCArr(cArr, ket, numQubits);
    // return trace_r(cArr, numQubits);

    double out = 0.0;
    for (auto const &el : ket)
    {
        out += el * el;
    }
    return out;
}

double fidel_ket(const std::vector<double> &in_ket, const std::vector<double> &corr_ket, const int numQubits)
{
    double bra_ket = 0.0;
    const int n = 1 << numQubits;
    for (int i = 0; i < n; i++)
    {
        bra_ket += in_ket[i] * corr_ket[i];
    }
    return bra_ket * bra_ket;
}

void disp_tr_fid(const std::vector<double> &in_ket, const int numQubits)
{
    std::cout << "fidelity = " << fidel_ket(in_ket, in_ket, numQubits) << std::endl;
    std::cout << "trace = " << trace_ket(in_ket, numQubits) << std::endl;
}

// Define the normalization factor
const double factor = 1.0 / std::sqrt(8.0);

// Function to print array as comma-separated values in curly braces
void printArray(const std::vector<double> &arr)
{
    std::cout << "{";
    for (size_t i = 0; i < arr.size(); ++i)
    {
        std::cout << arr[i];
        if (i < arr.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << "}" << std::endl;
}

int main()
{
    constexpr double theta = std::numbers::pi / 2;
    const double s = std::sin(0.5 * theta);
    const double c = std::cos(0.5 * theta);

    // Logical |0_L> codeword for the [[7, 1, 3]] Steane code
    std::vector<double> logical_ket = {c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0)};
    
    // ketX_r(logical_ket, 7, 2);
    // disp_tr_fid(logical_ket, 7);
    // ketY_r(logical_ket, 7, 3);
    // disp_tr_fid(logical_ket, 7);
    // ketZ_r(logical_ket, 7, 5);
    // disp_tr_fid(logical_ket, 7);

    std::cout << sizeof(double) << std::endl;
    
    // double num = 1 / std::numbers::sqrt2;
    // std::vector<double> ket_0 = {1.0, 0.0};
    // disp_tr_fid(ket_0, 1);
    // std::vector<double> ket_1 = {0.0, 1.0};
    // disp_tr_fid(ket_1, 1);
    // std::vector<double> ket_H1P = {num, num};
    // disp_tr_fid(ket_H1P, 1);
    // std::vector<double> ket_H1N = {num, -num};
    // disp_tr_fid(ket_H1N, 1);

    // std::vector<double> ket = {0.0, 1.0, 0.0, 0.0};
    // disp_tr_fid(ket, 2);
    // ketX_r(ket, 2, 0);
    // disp_tr_fid(ket, 2);

    // std::vector<double> ket2 = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // ketH_r(ket2, 3, 2);

    // std::vector<double> ket_01 = {0.0, 1.0, 0.0, 0.0};

    // std::vector<double> ket_p(8);
    // ketKron(ket_p, ket_H1P, ket_01, 1, 2);
    // disp_tr_fid(ket_p, 3);
    // ketCZ_r(ket_p, 3, 0, 2);
    // disp_tr_fid(ket_p, 3);

    // for (double value : ket_p) {
    //     std::cout << value << ", ";
    // }

    // std::cout << std::endl;

    // std::vector<double> out_ket;

    // applyAncillaKet(out_ket, 1, 3, ket_H1N, ket_p);

    // for (double value : out_ket) {
    //     std::cout << value << ", ";
    // }

    // std::vector<double> ket_p(4);
    // ketKron(ket_p, ket_H1P, ket_0, 1, 1);

    // ketCZ_r(ket_p, 2, 0, 1);

    // for (double value : ket_p) {
    //     std::cout << value << ", ";
    // }

    // double err = 0.01;

    // std::vector<double> eventProbabilities(13, err / 15);
    // eventProbabilities[0] = (1 - err) + 3 * err/15;

    // // std::cout << simulateMultipleEvents(eventProbabilities) << std::endl;

    // int i = 0, result = 1;

    // while(result == 1){
    //     result = simulateMultipleEvents(eventProbabilities);
    //     i++;
    // }

    // std::cout << result << " , " << i << std::endl;

    return 0;
}