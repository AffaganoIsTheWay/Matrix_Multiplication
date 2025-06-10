#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>

using namespace std;

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

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ofstream ResultFile;

    if(rank == 0){
        ResultFile.open("MPI_result.csv", ios_base::app);
    }

    int matrix_dimension[][3] = {{128, 128, 128},
                                 {256, 256, 256},
                                 {512, 512, 512},
                                 {1024, 1024, 1024},
                                 {2048, 2048, 2048},
                                 {4096, 4096, 4096},
                                 {2048, 128, 1024},
                                 {8192, 256, 128},
                                 {8192, 1024, 256}};

    for(int i = 0; i < 9; i++){
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
        int rows_per_proc = matrix_dimension[i][0] / size;
        int row_start = rank * rows_per_proc;
        int row_end = row_start + rows_per_proc;
        int* C_local = new int[rows_per_proc * matrix_dimension[i][2]];
        int* C = new int[matrix_dimension[i][0] * matrix_dimension[i][2]];

        for(int j = 0; j < 10; j++){

            if (rank == 0) {
                srand(time(NULL));
            
                for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][1]; ++k)
                    A[k] = rand() % 100 + 1;
            
                for (int k = 0; k < matrix_dimension[i][1] * matrix_dimension[i][2]; ++k)
                    B[k] = rand() % 100 + 1;
            
            }
        
            // Broadcast matrices A and B to all ranks
            MPI_Bcast(A, matrix_dimension[i][0] * matrix_dimension[i][1], MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(B, matrix_dimension[i][1] * matrix_dimension[i][2], MPI_INT, 0, MPI_COMM_WORLD);
        
        
            // Performing parallel multiplication
            MPI_Barrier(MPI_COMM_WORLD);
            double start_time = MPI_Wtime();

            MatrixMultiplicationMPI(A, B, C_local, matrix_dimension[i][0],
                                    matrix_dimension[i][1], matrix_dimension[i][2],
                                    row_start, row_end);

            // Ensure all process finished in order to gather the result
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Gather(C_local, rows_per_proc * matrix_dimension[i][2], MPI_INT,
                       C, rows_per_proc * matrix_dimension[i][2], MPI_INT,
                       0, MPI_COMM_WORLD);

            double end_time = MPI_Wtime();
            double duration = (end_time - start_time);

            if (rank == 0) {
                ResultFile << size << "," <<
                              matrix_dimension[i][0]  << "," <<
                              matrix_dimension[i][1] << "," <<
                              matrix_dimension[i][2] << "," <<
                              duration << endl;
            }
        }

        delete[] A;
        delete[] B;
        delete[] C_local;
        delete[] C;
    }

    if(rank == 0){
        ResultFile.close();
    }

    MPI_Finalize();
    return 0;
}
