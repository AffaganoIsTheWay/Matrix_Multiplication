#include <iostream>
#include <chrono>    // Time counter

using namespace std::chrono;

#define N 1000
//Define row and colums separate
//#define C1
//#define R1
//#define C2
//#define R2

float** matrixMultuply(float mat1[][N], float mat2[][N]){

    float** result = 0;
    result = new float*[N];

    for (int i = 0; i < N; i++) {
        result[i] = new float[N];
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
    float mat1[N][N], mat2[N][N];
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            mat1[i][j] = i + j;
            mat2[i][j] = i - j;
        }
    }

    //check if the multiplication is duable

    auto beg = high_resolution_clock::now();

    float** res = matrixMultuply(mat1, mat2);

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - beg);

    // Displaying the elapsed time
    std::cout << "Elapsed Time: " << duration.count() << " milliseconds" << std::endl;

    //bandwidth AND flops
    double data_transfer = 3.0 * ((N * N * N) * sizeof(int));
    double bandwidth = data_transfer / (duration.count() * 1e6);
    double flops = (N*N*N) / (duration.count() * 1e-3);

    std::cout << "Effective bandwith: " << bandwidth << " GB/s" << std::endl;
    std::cout << "FLOPS: " << flops << std::endl;
}
