#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <cmath>

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

void MatrixMultiplicationMPI(int* A, int* B, int* C,
                     int n, int m, int p,
                     int row_start, int row_end) {
    for (int i = row_start; i < row_end; ++i) {
        for (int j = 0; j < p; ++j) {
            int sum = 0;
            for (int k = 0; k < m; ++k) {
                sum += A[i * m + k] * B[k * p + j];
            }
            C[(i - row_start) * p + j] += sum;
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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    int P = atoi(argv[3]);

    if ((N % size != 0) || (M % size != 0) || (P % size != 0)) {
        if (rank == 0)
            cout << "The dimension must be divisible by the number of Threads\n";
        MPI_Finalize();
        return 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Inizialization
    int* A = new int[N * M];
    int* B = new int[M * P];

    // Partial C per rank
    int rows_per_proc = N / size;
    int row_start = rank * rows_per_proc;
    int row_end = row_start + rows_per_proc;
    int* C_local = new int[rows_per_proc * P];
    int* C_serial = nullptr;
    int* C_parallel = nullptr;

    if (rank == 0) {
        srand(time(NULL));
        for (int i = 0; i < N * M; ++i) A[i] = rand() % 100 + 1;
        for (int i = 0; i < M * P; ++i) B[i] = rand() % 100 + 1;

        // Performing Serial multiplication
        C_serial = new int[N * P];

        double start_serial = MPI_Wtime();

        MatrixMultiplication(A, B, C_serial, N, M, P);

        double end_serial = MPI_Wtime();
        double duration_serial = (end_serial - start_serial);
        double data_transferred_serial = 3.0 * (double)((double)N * (double)M * (double)P) * sizeof(int);
        double bandwidth_serial = data_transferred_serial / (duration_serial * 1e9);

        cout << "Time taken by serial: " << duration_serial << " seconds" << endl;
        cout << "Effective Serial Bandwidth: " << bandwidth_serial << " GB/s" << endl << endl;

        C_parallel = new int[N * P];
    }

    // Broadcast matrices A and B to all ranks
    MPI_Bcast(A, N * M, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, M * P, MPI_INT, 0, MPI_COMM_WORLD);

    // Performing parallel multiplication
    MPI_Barrier(MPI_COMM_WORLD);
    double start_parallel = MPI_Wtime();

    MatrixMultiplicationMPI(A, B, C_local, N, M, P, row_start, row_end);

    // Ensure all process finished in order to gather the result
    MPI_Barrier(MPI_COMM_WORLD);
    // Time needed just for the computation
    double end_computation = MPI_Wtime();

    MPI_Gather(C_local, rows_per_proc * P, MPI_INT,
               C_parallel, rows_per_proc * P, MPI_INT,
               0, MPI_COMM_WORLD);

    double end_parallel = MPI_Wtime();

    if (rank == 0) {

        double duration_parallel = (end_parallel - start_parallel);
        double duration_computation = (end_computation - start_parallel);
        double duration_gathering = (end_parallel - end_computation);
        double data_transferred_parallel = ((2.0 * (double)((double)N * (double)M * (double)P)) + (double)(N * P)) * sizeof(int);
        double bandwidth_parallel = data_transferred_parallel / (duration_parallel * 1e9);

        cout << "Time taken by parallel Rows Level: " << duration_parallel << " seconds\n"
             << " of which computation: " << duration_computation << endl
             << " and gathering results: " << duration_gathering << endl;
        cout << "Effective Parallel Bandwidth during computation: " << bandwidth_parallel << " GB/s" << endl << endl;

        cout << "Check multiplied Matrix:" << check_transpose(C_serial, C_parallel, N, P) << endl;

        delete[] C_serial;
        delete[] C_parallel;
    }

    delete[] A;
    delete[] B;
    delete[] C_local;

    MPI_Finalize();
    return 0;
}
