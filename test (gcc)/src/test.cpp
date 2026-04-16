// Just a file to test basic code

#include "helper.hpp"
#include <omp.h>

int main()
{
    constexpr double theta = std::numbers::pi / 2;
    const double s = std::sin(0.5 * theta);
    const double c = std::cos(0.5 * theta);

    // Logical |0_L> codeword for the [[7, 1, 3]] Steane code
    std::vector<double> logical_ket = {c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, c / std::sqrt(8.0), 0, 0, c / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, s / std::sqrt(8.0)};
    
    std::cout << logical_ket.max_size() << std::endl;

    int NUMTHREADS = omp_get_max_threads();
    std::cout << "Number of threads: " << NUMTHREADS << std::endl;
    // ketX_r(logical_ket, 7, 2);
    // disp_tr_fid(logical_ket, 7);
    // ketY_r(logical_ket, 7, 3);
    // disp_tr_fid(logical_ket, 7);
    // ketZ_r(logical_ket, 7, 5);
    // disp_tr_fid(logical_ket, 7);

    // std::cout << sizeof(double) << std::endl;
    
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