#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>

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

    if(rank == 0){
        ofstream ResultFile;
        ResultFile.open("MPI_result.csv", ios_base::app);
    }

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
        if ((matrix_dimension[i][0] % size != 0) || (matrix_dimension[i][1] % size != 0) || (matrix_dimension[i][2] % size != 0)) {
            if (rank == 0)
                cout << "The dimension must be divisible by the number of Threads\n";
            MPI_Finalize();
            return 1;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Inizialization
        int *A = new int[matrix_dimension[i][0] * matrix_dimension[i][1]];
        int *B = new int[matrix_dimension[i][1] * matrix_dimension[i][2]];
        // Partial C per rank
        int rows_per_proc = N / size;
        int row_start = rank * rows_per_proc;
        int row_end = row_start + rows_per_proc;
        int* C_local = new int[rows_per_proc * P];
        int* C_serial = nullptr;
        int* C_parallel = nullptr;

        for(int j = 0; j < 10; j++){
            if (rank == 0) {
                srand(time(NULL));

                for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][1]; ++k)
                    A[k] = rand() % 100 + 1;

                for (int k = 0; k < matrix_dimension[i][1] * matrix_dimension[i][2]; ++k)
                    B[k] = rand() % 100 + 1;

                // Performing Serial multiplication
                C_serial = new int[matrix_dimension[i][0] * matrix_dimension[i][2]];

                double start_serial = MPI_Wtime();

                MatrixMultiplication(A, B, C_serial, matrix_dimension[i][0], matrix_dimension[i][1], matrix_dimension[i][2]);

                double end_serial = MPI_Wtime();
                double duration_serial = (end_serial - start_serial);

                C_parallel = new int[matrix_dimension[i][0] * matrix_dimension[i][2]];
            }

            // Broadcast matrices A and B to all ranks
            MPI_Bcast(A, matrix_dimension[i][0] * matrix_dimension[i][1], MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(B, matrix_dimension[i][1] * matrix_dimension[i][2], MPI_INT, 0, MPI_COMM_WORLD);

            // Performing parallel multiplication
            MPI_Barrier(MPI_COMM_WORLD);
            double start_parallel = MPI_Wtime();

            MatrixMultiplicationMPI(A, B, C_local, matrix_dimension[i][0],
                                    matrix_dimension[i][1], matrix_dimension[i][2],
                                    row_start, row_end);

            // Ensure all process finished in order to gather the result
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Gather(C_local, rows_per_proc * P, MPI_INT,
                       C_parallel, rows_per_proc * P, MPI_INT,
                       0, MPI_COMM_WORLD);

            double end_parallel = MPI_Wtime();

            if (rank == 0) {
                ResultFile << omp_get_max_threads() << "," <<
                              matrix_dimension[i][0]  << "," <<
                              matrix_dimension[i][1] << "," <<
                              matrix_dimension[i][2] << "," <<
                              duration_serial << "," <<
                              duration_parallel << endl;
            }
        }

        delete[] A;
        delete[] B;
        delete[] C_local;
        if(rank == 0){
            delete[] C_serial;
            delete[] C_parallel;
        }
    }

    if(rank == 0){
        ResultFile.close();
    }

    MPI_Finalize();
    return 0;
}
