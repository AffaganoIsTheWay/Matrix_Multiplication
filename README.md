# IntroPARCO 2024 - BLAS GEMM

## Introduction

**BLAS GEMM (General Matrix-Matrix Multiplication)** is a fundamental routine in the Basic Linear Algebra Subprograms (BLAS) library used to perform high-performance matrix multiplication.

GEMM is widely used in scientific computing, machine learning, and graphics due to its efficiency and optimization on various hardware architectures.

For consistency throughout the document and to facilitate understanding, we define the two matrices as follows:

$$A = \begin{bmatrix} N \times M \end{bmatrix}$$  
$$B = \begin{bmatrix} M \times P \end{bmatrix}$$  
$$C = AB = \begin{bmatrix} N \times M \end{bmatrix} \cdot \begin{bmatrix} M \times P \end{bmatrix}$$

The performance of each method is measured by two key factors:

- **Execution Time**: Time taken to complete the matrix multiplication.
- **Memory Bandwidth**: Efficiency of memory usage during the multiplication.

---

To start clone the repository:

```bash
git clone https://github.com/AffaganoIsTheWay/Matrix_Multiplication.git
```

## Architecture

The tests in this repository were executed on the following machine architectures:

#### Personal Computers

- **Processor**: AMD Ryzen 5 5600G 3,90 GHz (6 cores, 12 threads)
- **Architecture**: Zen 3 (x86_64)
- **RAM**: 16 GB DDR4

<br>

- **Processor**: Intel(R) Core(R) i7-7500U CPU @ 2.70GHz (4 cores, 8 threads)
- **Architecture**: x86_64
- **RAM**: 16 GB DDR4

#### Cluster

- **Processor**: Intel(R) Xeon(R) Gold 6252N CPU @ 2.30GHz (96 cores, 96 threads)
- **Architecture**: x86_64
- **RAM**: 16 GB

---


## Sequential Implementation

In the sequential implementation, the matrix multiplication is performed by a single thread in a straightforward manner. The operation follows the standard matrix multiplication pattern, computing each element of the result matrix based on the dot product of corresponding rows and columns.

In this section, the **Execution Time** and **Memory Bandwidth** rapresents only the multiplication function.

To execute this measurament, first navigate to the directory:

```bash
cd Sequential\ Implementation/
```

Then run the bash script:

```bash
sh sequential_moltiplication.sh <number of iteration> <N> <M> <P>
```

Where:

- `<number of iterations>`: Number of times the test will be run.
- `<N>`, `<M>`, `<P>`: Dimensions of the matrices. For example, if `N=10`, `M=15`, `P=20`, then matrix *A* is 10×15 and *B* is 15×20.

Here a possible call:

```bash
sh sequential_moltiplication.sh 10 1024 1024 1024
```

---

## OpenMP

**OpenMP** is a popular API for parallel programming in C, C++, and Fortran. It enables easy parallelization of loops and sections of code by adding simple compiler directives.

In this section we will mesure the perfomances of different **OpenMP** implementation methods.

First, navigate to the directory

```bash
cd OpenMP/
```

You can run each implementation individually:

```bash
sh Implementation_script.sh <number of iterations> <number of threads> <N> <M> <P>
```

Or run all of them together:

```bash
sh all_compile_script.sh <number of iterations> <number of threads> <N> <M> <P>
```

Where:

- `<number of iterations>`: Number of test repetitions.
- `<number of threads>`: Number of threads to use for parallel execution.
- `<N>`, `<M>`, `<P>`: Matrix dimensions.

Here a possible call:

```bash
sh all_compile_script.sh 10 4 1024 1024 1024
```

The **Execution Time** and **Memory Bandwidth** are computed for both the sequential and the parallel methods and displayed.

---

## BLAS Library

Due to the widespread use of matrix multiplication in various fields, many libraries offer highly optimized implementations.

These optimizations include CPU-specific vector instructions (e.g., AVX, SSE), cache blocking to improve memory access patterns, multi-threading to utilize multiple cores, and algorithmic tuning to reduce computational overhead.

This project implemented the **OpenBlas** and the **Intel MKL** libraries

Navigate to the folder

```bash
cd BLAS\ library/
```

Run the implementation:

```bash
sh Implementation_script.sh <number of iterations> <number of threads> <N> <M> <P>
```

Where:

- `<number of iterations>`: Number of test repetitions.
- `<number of threads>`: Number of threads for parallel execution.
- `<N>`, `<M>`, `<P>`: Matrix dimensions.

Note: **Memory Bandwidth** is not measured in this section because the actual data access patterns are hidden within the library implementation.

### ⚠️ ATTENTION

Both libraries must be installed beforhand. In particular for the **Intel MKL** the `IntelBLAS_compile_script.sh` must be modified in order to match your `include` and `lib` directory paths

[Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-documentation.html)

[OpenBlas](http://www.openmathlib.org/OpenBLAS/)

---

## MPI

**MPI (Message Passing Interface)** is a standardized library used for writing parallel programs that run on multiple processors or computers. It allows these separate processes to communicate by sending and receiving messages, coordinating their work to solve large problems more efficiently. MPI is widely used in high-performance computing to build applications that can scale across clusters and supercomputers.

Navigato to the folder

```bash
cd MPI/
```

Run individual implementations:

```bash
sh Implementation_script.sh <number of iterations> <number of threads> <N> <M> <P>
```

Or run all implementations at once:

```bash
sh all_compile_script.sh <number of iterations> <number of threads> <N> <M> <P>
```

Parameters:

- `<number of iterations>`: Number of test repetitions.
- `<number of threads>`: Number of threads to run the simulation.
- `<N>`, `<M>`, `<P>`: Matrix dimensions.

Here a possible call:

```bash
sh all_compile_script.sh 10 4 1024 1024 1024
```

The **Execution Time** and **Memory Bandwidth** are computed for both the sequential and the parallel methods and displayed.

## Cluster

All main files perform multiplications between square matrices and tall/skinny matrices for multiple dimension sets (defined inside the files), each repeated 10 times. The results are saved in CSV files.

The `main_script.pbs` first of all ensure that all the required module are loaded and then call all the main files with different number of threads: `{1,2,4,8,16,32,64}`.

In order to run it, first, copy ```OMP_main.cpp```, ```MPI_ main.cpp```, ```BLASlib_main.cpp``` and ```main_script.pbs``` in the cluster. Open ```main_script.pbs``` and modify the working directory in to match yours. Then submit the job:

```bash
qsub main_script.pbs
```

A file named ```OMP_result.csv```, ```MPI_result.csv``` and ```BLASlib_result.csv``` will be created, in which all the result will be saved.

### ⚠️ ATTENTION

```main_script.pbs``` may require to be converted to unix format.

In order to do so, run the following command on the cluster:

```bash
dos2unix main_script.pbs
```
