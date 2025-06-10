#include <iostream>
#include <cstdlib>
#include <cblas.h>
#include <omp.h>
#include <fstream>

using namespace std;

int main(int argc, const char* argv[]) {
    
    ofstream ResultFile;
    ResultFile.open("BLASlib_result.csv", ios_base::app);

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
        // Inizialization matrix
        double *A = new double[matrix_dimension[i][0] * matrix_dimension[i][1]];
        double *B = new double[matrix_dimension[i][1] * matrix_dimension[i][2]];
        double *C = new double[matrix_dimension[i][0] * matrix_dimension[i][2]];

        for(int j = 0; j < 10; j++){
            srand(time(NULL));

            for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][1]; ++k)
                A[k] = (double)(rand() % 100 + 1);

            for (int k = 0; k < matrix_dimension[i][1] * matrix_dimension[i][2]; ++k)
                B[k] = (double)(rand() % 100 + 1);

            for (int k = 0; k < matrix_dimension[i][0] * matrix_dimension[i][2]; ++k)
                C[k] = 0.0;

            double start_time = omp_get_wtime();

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                matrix_dimension[i][1], matrix_dimension[i][0],
                matrix_dimension[i][2], 1.0, A, matrix_dimension[i][2],
                B, matrix_dimension[i][0], 0.0, C_parallel, matrix_dimension[i][0]);

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
        delete[] C_serial;
        delete[] C_parallel;
    }

    ResultFile.close();

    return 0;
}