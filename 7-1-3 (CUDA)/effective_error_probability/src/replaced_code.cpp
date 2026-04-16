    // Main Computation
    if (STRCOMPARE(mode, "flag"))
    {
#pragma omp parallel shared(logicalcArr, err, nLinksArray, ecInfo, rounds, logicalKet, errTeleFactor) private(etrans, j, correctedcArr, i, tempcArrCodeword) // Reference: https://stackoverflow.com/questions/2352895/how-to-ensure-a-dynamically-allocated-array-is-private-in-openmp
        {
#pragma omp for schedule(static, NUMPOINTS / NUMTHREADS)
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
        }
    }
    else
    {
#pragma omp parallel shared(logicalcArr, err, nLinksArray, ecInfo, rounds, logicalKet, errTeleFactor) private(etrans, j, correctedcArr, i, tempcArrCodeword) // Reference: https://stackoverflow.com/questions/2352895/how-to-ensure-a-dynamically-allocated-array-is-private-in-openmp
        {
#pragma omp for schedule(static, NUMPOINTS / NUMTHREADS)
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
                    correctedcArr = QubitCode::OneRoundSequential(tempcArrCodeword,
                                                                      err / errReencodingFactor,
                                                                      errTeleFactor * err / errReencodingFactor,
                                                                      QubitCode::correctionMode::oneLoss,
                                                                      0);
                    // NOTE: Is this "0" above really appropriate without Reset?
                    ecInfo[j].grandErrors[i] = (1 - FastPauli::fidel_r(logicalKet, correctedcArr, NUMQUBITS));
                    ecInfo[j].grandTraces[i] = FastPauli::trace_r(correctedcArr, NUMQUBITS);

                    tempcArrCodeword = correctedcArr;
                }
            }
        }
    }