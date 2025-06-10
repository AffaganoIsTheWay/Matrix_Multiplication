#include <iostream>
#include <cstdlib>
#include <cblas.h>
#include <omp.h>
#include <fstream>

using namespace std;

void MatrixMultiplication(double* A, double* B, double* C, int n, int m, int p) {
    // Multiply matrices
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < m; ++k) {
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}

bool check_transpose(double* C_serial,double* C_parallel, int N, int P){
    for (int i = 0; i < N * P; ++i) {
        if(C_serial[i] != C_parallel[i]){
            return false;
        }
    }

    return true;
}

int main(int argc, const char* argv[]) {
    
    ofstream ResultFile;
    ResultFile.open("BLASlib_result.csv", ios_base::app);

    int matrix_dimension[][3] = {{128, 128, 128},
                                 {256, 256, 256},
                                 {512, 512, 512},
                                 {1024, 1024, 1024},
                                 {2048, 2048, 2048},
                                 {4096, 4096, 4096},
                                 {8192, 8192, 8192},
                                 {2048, 128, 1024},
                                 {8192, 256, 128},
                                 {8192, 1024, 256}};

    for(int i = 0; i < 10; i++){
        // Inizialization matrix
        double *A = new double[matrix_dimension[i][0] * matrix_dimension[i][1]];
        double *B = new double[matrix_dimension[i][1] * matrix_dimension[i][2]];
        double *C_serial = new double[matrix_dimension[i][0] * matrix_dimension[i][2]];
        double *C_parallel = new double[matrix_dimension[i][0] * matrix_dimension[i][2]];

        for(int j = 0; j < 10; j++){
            srand(time(NULL));

            for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][1]; ++k)
                A[k] = (double)(rand() % 100 + 1);

            for (int k = 0; k < matrix_dimension[i][1] * matrix_dimension[i][2]; ++k)
                B[k] = (double)(rand() % 100 + 1);

            for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][2]; ++k){
                C_serial[k] = 0.0;
                C_parallel[k] = 0.0;
            }

            double start_serial = omp_get_wtime();

            MatrixMultiplication(A, B, C_serial, matrix_dimension[i][0], matrix_dimension[i][1], matrix_dimension[i][2]);

            double end_serial = omp_get_wtime();
            double duration_serial = (end_serial - start_serial);

            double start_parallel = omp_get_wtime();

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                matrix_dimension[i][1], matrix_dimension[i][0],
                matrix_dimension[i][2], 1.0, A, matrix_dimension[i][2],
                B, matrix_dimension[i][0], 0.0, C_parallel, matrix_dimension[i][0]);

            double end_parallel = omp_get_wtime();
            double duration_parallel = (end_parallel - start_parallel);

            ResultFile << openblas_get_max_threads() << "," <<
                          matrix_dimension[i][0]  << "," <<
                          matrix_dimension[i][1] << "," <<
                          matrix_dimension[i][2] << "," <<
                          duration_serial << "," <<
                          duration_parallel << endl;
        }

        // Cleanup memory
        delete[] A;
        delete[] B;
        delete[] C_serial;
        delete[] C_parallel;
    }

    ResultFile.close();

    return 0;
}