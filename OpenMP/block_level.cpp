#include <iostream>
#include <omp.h>
#include <ctime>
#include <cstdlib>

using namespace std;

void MatrixMultiplication(int* A, int* B, int* C, int n, int m, int p) {
    // Multiply matrices
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < m; ++k) {
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}

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

bool check_transpose(int* C_serial,int* C_parallel, int N, int P){
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
    int* A = new int[N * M];
    int* B = new int[M * P];
    int* C_serial = new int[N * P];
    int* C_parallel = new int[N * P];
    int blockSize = 64;

    srand(time(NULL));

    for (int i = 0; i < N * M; ++i) A[i] = rand() % 100 + 1;
    for (int i = 0; i < M * P; ++i) B[i] = rand() % 100 + 1;

    for (int i = 0; i < N * P; ++i) {
        C_serial[i] = 0;
        C_parallel[i] = 0;
    }

    /* double start_serial = omp_get_wtime();

    MatrixMultiplication(A, B, C_serial, N, M, P);

    double end_serial = omp_get_wtime();
    double duration_serial = (end_serial - start_serial);
    double data_transferred_serial = 3.0 * (double)((double)N * (double)M * (double)P) * sizeof(int);
    double bandwidth_serial = data_transferred_serial / (duration_serial * 1e9);

    cout << "Time taken by serial: " << duration_serial << " seconds" << endl;
    cout << "Effective Serial Bandwidth: " << bandwidth_serial << " GB/s" << endl << endl; */

    double start_parallel = omp_get_wtime();

    MatrixMultiplicationOMP(A, B, C_parallel, N, M, P, blockSize);

    double end_parallel = omp_get_wtime();
    double duration_parallel = (end_parallel - start_parallel);
    double data_transferred_parallel = ((2.0 * (double)((double)N * (double)M * (double)P)) + (double)(N * P)) * sizeof(int);
    double bandwidth_parallel = data_transferred_parallel / (duration_parallel * 1e9);

    cout << "Time taken by parallel Block Level: " << duration_parallel << " seconds" << endl;
    cout << "Effective Parallel Bandwidth: " << bandwidth_parallel << " GB/s" << endl << endl;

    cout << "Check multiplied Matrix:" << check_transpose(C_serial, C_parallel, N, P) << endl;

    // Cleanup memory
    delete[] A;
    delete[] B;
    delete[] C_serial;
    delete[] C_parallel;

    return 0;
}