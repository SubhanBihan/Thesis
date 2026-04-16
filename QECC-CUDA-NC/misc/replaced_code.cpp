std::vector<cuda_global> X_d;
std::vector<cuda_global> Y_d;
std::vector<cuda_global> Z_d;

    Dynamic_t tmp_XYZ_h, I_i{1};
    for (int i = 0; i < NUMQUBITS; i++)
    {
        Dynamic_t X_i{1}, Y_i{1}, Z_i{1};

        I_i = kroneckerProduct(I_i, I_h).eval();
        for (int j = 0; j < NUMQUBITS; j++) {
            X_i = kroneckerProduct(X_i, (i == j) ? X_h : I_h).eval();
            Y_i = kroneckerProduct(Y_i, (i == j) ? Y_h : I_h).eval();
            Z_i = kroneckerProduct(Z_i, (i == j) ? Z_h : I_h).eval();
        }

        X_d.push_back(cuda_global(X_i));
        Y_d.push_back(cuda_global(Y_i));
        Z_d.push_back(cuda_global(Z_i));
    }

    cuda_global I_d(I_i);



CHECK_CUBLAS(cublasSgemm(handle1, CUBLAS_OP_N, CUBLAS_OP_N, DIM, DIM, DIM, &alpha, X_d[k].get(), DIM, tempDM_d.get(), DIM, &beta,  tempX1_d.get(), DIM));
CHECK_CUBLAS(cublasSgemm(handle2, CUBLAS_OP_N, CUBLAS_OP_N, DIM, DIM, DIM, &alpha, Y_d[k].get(), DIM, tempDM_d.get(), DIM, &beta,  tempY1_d.get(), DIM));
CHECK_CUBLAS(cublasSgemm(handle3, CUBLAS_OP_N, CUBLAS_OP_N, DIM, DIM, DIM, &alpha, Z_d[k].get(), DIM, tempDM_d.get(), DIM, &beta,  tempZ1_d.get(), DIM));

CHECK_CUBLAS(cublasSgemm(handle1, CUBLAS_OP_N, CUBLAS_OP_N, DIM, DIM, DIM, &alpha, tempX1_d.get(), DIM, X_d[k].get(), DIM, &beta,  tempX2_d.get(), DIM));
CHECK_CUBLAS(cublasSgemm(handle2, CUBLAS_OP_N, CUBLAS_OP_N, DIM, DIM, DIM, &alpha, tempY1_d.get(), DIM, Y_d[k].get(), DIM, &beta,  tempY2_d.get(), DIM));
CHECK_CUBLAS(cublasSgemm(handle3, CUBLAS_OP_N, CUBLAS_OP_N, DIM, DIM, DIM, &alpha, tempZ1_d.get(), DIM, Z_d[k].get(), DIM, &beta,  tempZ2_d.get(), DIM));
