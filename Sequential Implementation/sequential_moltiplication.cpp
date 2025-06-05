#include <iostream>
#include <omp.h>
#include <ctime>
#include <cstdlib>

using namespace std;

void multiplyMatrices(int* A, int* B, int* C, int n, int m, int p) {
    // Multiply matrices
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < m; ++k) {
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}

int main(int argc, const char* argv[]) {
    // Dimensions: A is n x m, B is m x p
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    int P = atoi(argv[3]);

    // Inizialization matrix
    int* A = new int[N * M];
    int* B = new int[M * P];
    int* C = new int[N * P];

    srand(time(NULL));

    for (int i = 0; i < N; ++i){
        for (int j = 0; j < M; ++j) {
            A[i * M + j] = rand() % 100 + 1;
        }
    }

    for (int i = 0; i < M; ++i){
        for (int j = 0; j < P; ++j) {
            B[i * P + j] = rand() % 100 + 1;
        }
    }

    for (int i = 0; i < N * P; ++i) {
        C[i] = 0;
    }

    double start = omp_get_wtime();

    multiplyMatrices(A, B, C, N, M, P);

    double end = omp_get_wtime();
    double duration = (end - start);
    double data_transferred = 3 * (N * M * P) * sizeof(int);
    double bandwidth = data_transferred / (duration * 1e9);

    std::cout << "Time taken by serial: " << duration << " seconds" << endl;
    std::cout << "Effective Serial Bandwidth: " << bandwidth << " GB/s" << endl << endl;

    // Cleanup memory
    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}