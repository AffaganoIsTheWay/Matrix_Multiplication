#include <iostream>
#include <cstdlib>
#include <cblas.h>
#include <omp.h>

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
    // Dimensions: A is n x m, B is m x p
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    int P = atoi(argv[3]);

    // Inizialization matrix
    double* A = new double[N * M];
    double* B = new double[M * P];
    double* C_serial = new double[N * P];
    double* C_parallel = new double[N * P];

    srand(time(NULL));

    for (int i = 0; i < N * M; ++i) A[i] = (double)(rand() % 100 + 1);
    for (int i = 0; i < M * P; ++i) B[i] = (double)(rand() % 100 + 1);

    for (int i = 0; i < N * P; ++i) {
        C_serial[i] = 0.0;
        C_parallel[i] = 0.0;
    }

    double start_serial = omp_get_wtime();

    MatrixMultiplication(A, B, C_serial, N, M, P);

    double end_serial = omp_get_wtime();
    double duration_serial = (end_serial - start_serial);
    double data_transferred_serial = 3.0 * (double)((double)N * (double)M * (double)P) * sizeof(int);
    double bandwidth_serial = data_transferred_serial / (duration_serial * 1e9);

    cout << "Time taken by serial: " << duration_serial << " seconds" << endl;
    cout << "Effective Serial Bandwidth: " << bandwidth_serial << " GB/s" << endl << endl;

    double start_parallel = omp_get_wtime();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                N, P, M, 1.0, A, M, B, P, 0.0, C_parallel, P);
    
                double end_parallel = omp_get_wtime();
    double duration_parallel = (end_parallel - start_parallel);

    cout << "Time taken by OpenBLAS: " << duration_parallel << " seconds" << endl;

    cout << "Check multiplied Matrix:" << check_transpose(C_serial, C_parallel, N, P) << endl;

    // Cleanup memory
    delete[] A;
    delete[] B;
    delete[] C_serial;
    delete[] C_parallel;

    return 0;
}