// Just a file to test basic code

#include <iostream>
#include <filesystem>
#include <string>
#include <format>
#include <vector>

std::string findFileWithPartName(const std::string& directoryPath, const std::string& partOfFileName) {
    namespace fs = std::filesystem;

    // Iterate through the directory
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        // Check if the entry is a file
        if (entry.is_regular_file()) {
            // Get the filename as a string
            std::string filename = entry.path().filename().string();

            // Check if the filename contains the partOfFileName
            if (filename.find(partOfFileName) != std::string::npos) {
                // Return the full path of the file
                return entry.path().string();
            }
        }
    }

    // Return an empty string if no file is found
    return "";
}

void apply_permutation(std::vector<double>& v, int *indices, const int vSize) {  // Reference: https://stackoverflow.com/questions/67751784/how-to-do-in-place-sorting-a-list-according-to-a-given-index-in-c.
    for (int i = 0; i < vSize; i++) {
        auto current = i;
        while (i != indices[current]) {
            auto next = indices[current];
            std::swap(v[current], v[next]);
            indices[current] = current;
            current = next;
        }
        indices[current] = current;
    }
}

void ketX_r(std::vector<double>& ket, const int numQubits, const int targetQubit) {
    int i;
    const int n = 1 << numQubits;  // Matrix dimension.
    const int tmp0 = (n >> (targetQubit + 1));
    int *ketIndex = new int[n];
    for (i = 0; i < n; i++) {
        ketIndex[i] = i ^ tmp0;
    }
    apply_permutation(ket, ketIndex, n);
    delete[] ketIndex;
}

int main() {
    // std::string path = "../effective_error_probability/binary_data/flag_err1e-03_reencodingFactor3_teleFactor3";
    // std::string partOfFileName = "ecInfo0";

    // std::string foundFile = findFileWithPartName(path, partOfFileName);
    // if (!foundFile.empty()) {
    //     std::cout << "File found: " << foundFile << std::endl;
    // } else {
    //     std::cout << "No file found with the specified part of the name." << std::endl;
    // }
    
    // int i = 2;
    // std::cout << std::format("\x1b[0;92mecInfo.traces[{}] = {", i);

    std::vector<double> ket = {1.0, 0.0, 0.0, 0.0};
    ketX_r(ket, 2, 0);

    for (double value : ket) {
        std::cout << value << ", ";
    }

    return 0;
}
