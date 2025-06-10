#include <iostream>
#include <omp.h>
#include <ctime>
#include <cstdlib>
#include <fstream>

using namespace std;

void MatrixMultiplicationOMP(int* A, int* B, int* C, int n, int m, int p, int blockSize) {
    // Multiply matrices block parallelization
    #pragma omp parallel for collapse(2)
    for (int ii = 0; ii < n; ii += blockSize) {
        for (int jj = 0; jj < p; jj += blockSize) {
            for (int kk = 0; kk < m; kk += blockSize) {
                for (int i = ii; i < std::min(ii + blockSize, n); ++i) {
                    for (int j = jj; j < std::min(jj + blockSize, p); ++j) {
                        int sum = 0;
                        for (int k = kk; k < std::min(kk + blockSize, m); ++k) {
                            sum += A[i * m + k] * B[k * p + j];
                        }
                        #pragma omp atomic
                        C[i * p + j] += sum;
                    }
                }
            }
        }
    }
}

int main(int argc, const char *argv[]){

    ofstream ResultFile;
    ResultFile.open("OMP_result.csv", ios_base::app);

    int matrix_dimension[][3] = {{128, 128, 128},
                                 {256, 256, 256},
                                 {512, 512, 512},
                                 {1024, 1024, 1024},
                                 {2048, 2048, 2048},
                                 {4096, 4096, 4096},
                                 {2048, 128, 1024},
                                 {8192, 256, 128},
                                 {8192, 1024, 256}};

    for (int i = 0; i < 9; i++){
        // Inizialization matrix
        int *A = new int[matrix_dimension[i][0] * matrix_dimension[i][1]];
        int *B = new int[matrix_dimension[i][1] * matrix_dimension[i][2]];
        int *C = new int[matrix_dimension[i][0] * matrix_dimension[i][2]];
        int blockSize = 64;

        for (int j = 0; j < 10; j++){
            srand(time(NULL));

            for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][1]; ++k)
                A[k] = rand() % 100 + 1;

            for (int k = 0; k < matrix_dimension[i][1] * matrix_dimension[i][2]; ++k)
                B[k] = rand() % 100 + 1;

            for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][2]; ++k)
                C[k] = 0;

            double start_time = omp_get_wtime();

            MatrixMultiplicationOMP(A, B, C, matrix_dimension[i][0], matrix_dimension[i][1], matrix_dimension[i][2], blockSize);

            double end_time = omp_get_wtime();
            double duration = (end_time - start_time);

            ResultFile << omp_get_max_threads() << "," <<
                          matrix_dimension[i][0]  << "," <<
                          matrix_dimension[i][1] << "," <<
                          matrix_dimension[i][2] << "," <<
                          duration << endl;
        }

        // Cleanup memory
        delete[] A;
        delete[] B;
        delete[] C;
    }

    ResultFile.close();

    return 0;
}