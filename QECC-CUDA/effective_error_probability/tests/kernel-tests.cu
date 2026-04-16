#include "AuxMath.cuh"

#include "prereq.cuh"

void kernel_tests () {
    Matrix X(2, 2, {0, 1, 1, 0}), Y(2, 2, {0, 1, -1, 0}), Z(2, 2, {1, 0, 0, -1});
    unsigned m = 2, n = 1;
    Matrix ket(m, n), bra(n, m), DM(m, m);
    Matrix XpX_CPU(m ,m), YpY_CPU(m ,m), ZpZ_CPU(m ,m);
    Matrix XpX_cGPU(m ,m), YpY_cGPU(m ,m), ZpZ_cGPU(m ,m);
    cuda_global XpX_GPU(m * m), YpY_GPU(m * m), ZpZ_GPU(m * m);

    std::vector<std::vector<cudaStream_t>> streams_sets(2, std::vector<cudaStream_t>(NUM_STREAMS));

    for (auto &streams : streams_sets) {
        for (auto &stream : streams) {
            CHECK_CUDA(cudaStreamCreate(&stream));
        }
    }

    for(float theta : AuxMath::linspace(0, 2*pi, 10000)) {
        ket.mat = {cos(theta), sin(theta)};
        bra = AuxMath::bra(ket);
        DM = AuxMath::kronecker_col_major(ket, bra);

        XpX_CPU = X * DM * X;
        YpY_CPU = Y * DM * Y;
        ZpZ_CPU = Z * DM * Z;

        cuda_global DM_d(DM.mat);

        X_dm<<<gridSize, blockSize, 0, streams_sets[0][1]>>>(XpX_GPU.get(), DM_d.get(), 0, DIM);

        Y_dm<<<gridSize, blockSize, 0, streams_sets[0][2]>>>(YpY_GPU.get(), DM_d.get(), 0, DIM);

        Z_dm<<<gridSize, blockSize, 0, streams_sets[0][3]>>>(ZpZ_GPU.get(), DM_d.get(), 0, DIM);

        SyncStreams(streams_sets[0]);

        XpX_GPU.copy_to_host(XpX_cGPU.mat);
        YpY_GPU.copy_to_host(YpY_cGPU.mat);
        ZpZ_GPU.copy_to_host(ZpZ_cGPU.mat);

        if (!AuxMath::is_approx_equal(XpX_CPU, XpX_cGPU, 1e-5) ||
            !AuxMath::is_approx_equal(YpY_CPU, YpY_cGPU, 1e-5) ||
            !AuxMath::is_approx_equal(ZpZ_CPU, ZpZ_cGPU, 1e-5)) {
            std::cerr << "Test failed at theta = " << theta << "\n";
            std::cerr << "XpX_CPU:\n";
            print1Dvec(XpX_CPU); 
            std::cerr << "XpX_GPU:\n";
            print1Dvec(XpX_cGPU);
            std::cerr << "YpY_CPU:\n";
            print1Dvec(YpY_CPU); 
            std::cerr << "YpY_GPU:\n";
            print1Dvec(YpY_cGPU);
            std::cerr << "ZpZ_CPU:\n";
            print1Dvec(ZpZ_CPU); 
            std::cerr << "ZpZ_GPU:\n";
            print1Dvec(ZpZ_cGPU);
            
            std::exit(EXIT_FAILURE);
        }

    }

}