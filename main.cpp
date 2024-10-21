#include <iostream>
#include <ctime>
#include <cstdlib>

#define N 100
//Define row and colums separate
//#define C1
//#define R1
//#define C2
//#define R2

float** matrixMultuply(float** mat1, float** mat2){

    float** result = static_cast<float**>(malloc(N * sizeof(float)));

    for (int i = 0; i < N; i++) {
        result[i] = static_cast<float*>(malloc(N * sizeof(float)));
        for (int j = 0; j < N; j++) {
            result[i][j] = 0;

            for (int k = 0; k < N; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}

int main(){

    float** mat1 = static_cast<float**>(malloc(N * sizeof(float)));
    float** mat2 = static_cast<float**>(malloc(N * sizeof(float)));

    for(int i = 0; i < N; i++){
        mat1[i] = static_cast<float*>(malloc(N * sizeof(float)));
        mat2[i] = static_cast<float*>(malloc(N * sizeof(float)));
        for(int j = 0; j < N; j++){
            mat1[i][j] = i + j;
            mat2[i][j] = i - j;
        }
    }

    //check if the multiplication is duable

    clock_t start = clock();
    
    float** res = matrixMultuply(mat1, mat2);

    clock_t end = clock();
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;  // Time in seconds

    // Displaying the elapsed time
    std::cout << "Elapsed Time: " << duration << " seconds" << std::endl;

    //bandwidth AND flops
    double data_transfer = 3.0 * ((N * N * N) * sizeof(int));
    double bandwidth = data_transfer / (duration * 1e8);
    double flops = (N*N*N) / duration;

    std::cout << "Effective bandwith: " << bandwidth << " GB/s" << std::endl;
    std::cout << "FLOPS: " << flops << std::endl;

    free(mat1);
    free(mat2);
    free(res);

    return 0;
}
