# IntroPARCO 2024 - BLAS GEMM

## Introduction

**BLAS GEMM (General Matrix-Matrix Multiplication)** is a fundamental routine in the Basic Linear Algebra Subprograms (BLAS) library used to perform high-performance matrix multiplication.

GEMM is widely used in scientific computing, machine learning, and graphics due to its efficiency and optimization on various hardware architectures.

#

For cosistency inside the document and for a easier understanding we now define the two matrices.

$$A = \begin{bmatrix} N \times M \end{bmatrix}$$

$$B = \begin{bmatrix} M \\ \times \\ P \end{bmatrix}$$

And the operation performed.

$$C = AB = \begin{bmatrix} N \times M \end{bmatrix} \cdot \begin{bmatrix} M \\ \times \\ P \end{bmatrix}$$

#

The performances of every method are mesured by this two key factors:

- **Execution Time**: Time taken to complete the matrix multiplication.
- **Memory Bandwidth**: The efficiency of memory usage during the multiplication.

#

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
- **RAM**: 8 GB

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

- **number of itaration** is how many times the test will be performed.
- **N, M and P** are the sizes of the matrices (if 10, 15 and 20 the matrices will be *A*=10x15 and *B**A*=10x15 and *B*=15x20).

Here a possible call:

```bash
sh sequential_moltiplication.sh 10 1024 1024 1024
```

## OpenMP

OpenMP is a popular API for parallel programming in C, C++, and Fortran. It enables easy parallelization of loops and sections of code by adding simple compiler directives.

In this section different OpenMP implementation are compared in order to find the best solution.

#

In this section we will mesure the perfomances of different **OpenMP** implementation methods.

The **Execution Time** and **Memory Bandwidth** are computed for both the sequential and the parallel methods and displayed.

First, navigate to the directory

```bash
cd OpenMP/
```

Here we can run each of the implementation on its own,

```bash
sh Implementation_script.sh <number of iteration> <number of threads> <N> <M> <P>
```

Or we can all of them together with:

```bash
sh all_compile_script.sh <number of iteration> <number of threads> <N> <M> <P>
```
- **number of itaration** is how many times the test will be performed.
- **number of threads** is with how many threads you want to run the simulation with.
- **N, M and P** are the sizes of the matrices (if 10, 15 and 20 the matrices will be *A*=10x15 and *B**A*=10x15 and *B*=15x20).

Here a possible call:

```bash
sh all_compile_script.sh 10 4 1024 1024 1024
```

## BLAS Library

Due to the vast usage of the matrix multiplication in varius filds, there exist some library that contains an already hight optimized implementation.

These optimizations include CPU-specific vector instructions (e.g., AVX, SSE), cache blocking to improve memory access patterns, multi-threading to utilize multiple cores, and algorithmic tuning to reduce computational overhead.

For this project I implemented the **OpenBlas** and the **Intel MKL** libraries

#

Navigate to the folder

```bash
cd BLAS\ library/
```

Here we can run the implementations:

```bash
sh Implementation_script.sh <number of iteration> <number of threads> <N> <M> <P>
```

- **number of itaration** is how many times the test will be performed.
- **number of threads** is with how many threads you want to run the simulation with.
- **N, M and P** are the sizes of the matrices (if 10, 15 and 20 the matrices will be *A*=10x15 and *B**A*=10x15 and *B*=15x20).

### ATTENTION

Both this libraries require to be installed beforhand. In particular for the **Intel MKL** the ```IntelBLAS_compile_script.sh``` must be modified in order to match your *include* and *lib* directory paths

[Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-documentation.html)

[OpenBlas](http://www.openmathlib.org/OpenBLAS/)

## MPI

**MPI (Message Passing Interface)** is a standardized library used for writing parallel programs that run on multiple processors or computers. It allows these separate processes to communicate by sending and receiving messages, coordinating their work to solve large problems more efficiently. MPI is widely used in high-performance computing to build applications that can scale across clusters and supercomputers.

#

Navigato to the folder

```bash
cd MPI/
```

Here we can run each of the implementation on its own,

```bash
sh Implementation_script.sh <number of iteration> <number of threads> <N> <M> <P>
```

Or we can all of them together with:

```bash
sh all_compile_script.sh <number of iteration> <number of threads> <N> <M> <P>
```
- **number of itaration** is how many times the test will be performed.
- **number of threads** is with how many threads you want to run the simulation with.
- **N, M and P** are the sizes of the matrices (if 10, 15 and 20 the matrices will be *A*=10x15 and *B**A*=10x15 and *B*=15x20).

Here a possible call:

```bash
sh all_compile_script.sh 10 4 1024 1024 1024
```