#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

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

bool check_transpose(int* C_serial,int* C_parallel, int N, int P){
    for (int i = 0; i < N * P; ++i) {
        if(C_serial[i] != C_parallel[i]){
            return false;
        }
    }

    return true;
}

void MatrixMultiplicationMPI(int* A, int* B, int* C,
                     int n, int m, int p,
                     int blockN, int blockP, int blockM) {
    for (int i = 0; i < blockN; ++i) {
        for (int j = 0; j < blockP; ++j) {
            int sum = 0;
            for (int k = 0; k < blockM; ++k) {
                sum += A[i * blockM + k] * B[k * blockP + j];
            }
            C[i * blockP + j] = sum;
        }
    }
}


int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    int P = atoi(argv[3]);

    // Determine process grid
    int dims = std::sqrt(size);
    if (dims * dims != size) {
        if (rank == 0) cout << "Error: number of processes must be a perfect square." << endl;
        MPI_Finalize();
        return 1;
    }

    if (N % dims != 0 || M % dims != 0 || P % dims != 0) {
        if (rank == 0)
            cout << "Error: N, M, and P must be divisible by sqrt(num_processes)." << endl;
        MPI_Finalize();
        return 1;
    }

    int blockN = N / dims;
    int blockM = M / dims;
    int blockP = P / dims;

    int row_rank = rank / dims;
    int col_rank = rank % dims;

    // Allocate full matrices on root
    int* A_full = nullptr;
    int* B_full = nullptr;
    int* C_serial = nullptr;
    int* C_parallel = nullptr;

    if (rank == 0) {
        A_full = new int[N * M];
        B_full = new int[M * P];
        C_serial = new int[N * P];
        C_parallel = new int[N * P];

        srand(time(NULL));
        for (int i = 0; i < N * M; ++i) A_full[i] = rand() % 100 + 1;
        for (int i = 0; i < M * P; ++i) B_full[i] = rand() % 100 + 1;

        // Performing Serial multiplication
        C_serial = new int[N * P];

        double start_serial = MPI_Wtime();

        MatrixMultiplication(A_full, B_full, C_serial, N, M, P);
        
        double end_serial = MPI_Wtime();
        double duration_serial = (end_serial - start_serial);
        double data_transferred_serial = 3.0 * (double)((double)N * (double)M * (double)P) * sizeof(int);
        double bandwidth_serial = data_transferred_serial / (duration_serial * 1e9);

        cout << "Time taken by serial: " << duration_serial << " seconds" << endl;
        cout << "Effective Serial Bandwidth: " << bandwidth_serial << " GB/s" << endl << endl;
    }

    // Allocate local blocks
    int* A_block = new int[blockN * blockM];
    int* B_block = new int[blockM * blockP];
    int* C_block = new int[blockN * blockP]();
    
    // Scatter A blocks by row
    int* A_row_block = new int[blockN * M];
    if (col_rank == 0) {
        // Root distributes row blocks of A
        if (rank == 0) {
            for (int i = 0; i < dims; ++i) {
                for (int j = 0; j < blockN; ++j) {
                    int global_row = i * blockN + j;
                    for (int k = 0; k < M; ++k) {
                        A_row_block[j * M + k] = A_full[global_row * M + k];
                    }
                }
                if (i == 0)
                    std::copy(A_row_block, A_row_block + blockN * M, A_row_block);
                else
                    MPI_Send(A_row_block, blockN * M, MPI_INT, i * dims, 0, MPI_COMM_WORLD);
            }
        } else if (rank % dims == 0) {
            MPI_Recv(A_row_block, blockN * M, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // Each process gets its A_block from A_row_block
    MPI_Bcast(A_row_block, blockN * M, MPI_INT, row_rank * dims, MPI_COMM_WORLD);
    for (int i = 0; i < blockN; ++i)
        for (int j = 0; j < blockM; ++j)
            A_block[i * blockM + j] = A_row_block[i * M + col_rank * blockM + j];

    // Scatter B blocks by column
    int* B_col_block = new int[M * blockP];
    if (row_rank == 0) {
        if (rank == 0) {
            for (int j = 0; j < dims; ++j) {
                for (int i = 0; i < M; ++i) {
                    for (int k = 0; k < blockP; ++k) {
                        B_col_block[i * blockP + k] = B_full[i * P + j * blockP + k];
                    }
                }
                if (j == 0)
                    std::copy(B_col_block, B_col_block + M * blockP, B_col_block);
                else
                    MPI_Send(B_col_block, M * blockP, MPI_INT, j, 1, MPI_COMM_WORLD);
            }
        } else if (rank < dims) {
            MPI_Recv(B_col_block, M * blockP, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    MPI_Bcast(B_col_block, M * blockP, MPI_INT, col_rank, MPI_COMM_WORLD);
    for (int i = 0; i < blockM; ++i)
        for (int j = 0; j < blockP; ++j)
            B_block[i * blockP + j] = B_col_block[(row_rank * blockM + i) * blockP + j];

    MPI_Barrier(MPI_COMM_WORLD);
    double start_parallel = MPI_Wtime();

    MatrixMultiplicationMPI(A_block, B_block, C_block, N, M, P, blockN, blockP, blockM);

    // Ensure all process finished in order to gather the result
    MPI_Barrier(MPI_COMM_WORLD);
    // Time needed just for the computation
    double end_computation = MPI_Wtime();

    // Gather C blocks to rank 0
    MPI_Gather(C_block, blockN * blockP, MPI_INT,
               C_parallel, blockN * blockP, MPI_INT,
               0, MPI_COMM_WORLD);

    double end_parallel = MPI_Wtime();

    if (rank == 0) {

        double duration_parallel = (end_parallel - start_parallel);
        double duration_computation = (end_computation - start_parallel);
        double duration_gathering = (end_parallel - end_computation);
        double data_transferred_parallel = ((2.0 * (double)((double)N * (double)M * (double)P)) + (double)(N * P)) * sizeof(int);
        double bandwidth_parallel = data_transferred_parallel / (duration_parallel * 1e9);

        cout << "Time taken by parallel Block Level: " << duration_parallel << " seconds\n"
             << " of which computation: " << duration_computation << endl
             << " and gathering results: " << duration_gathering << endl;
        cout << "Effective Parallel Bandwidth during computation: " << bandwidth_parallel << " GB/s" << endl << endl;

        cout << "Check multiplied Matrix:" << check_transpose(C_serial, C_parallel, N, P) << endl;
    }

    delete[] A_block;
    delete[] B_block;
    delete[] C_block;
    delete[] A_row_block;
    delete[] B_col_block;
    if (rank == 0) {
        delete[] A_full;
        delete[] B_full;
        delete[] C_serial;
        delete[] C_parallel;
    }

    MPI_Finalize();
    return 0;
}