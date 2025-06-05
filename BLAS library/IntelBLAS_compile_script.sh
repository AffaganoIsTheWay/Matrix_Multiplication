#!/bin/bash

# Check argument
if [ $# -lt 5 ]; then
  echo "Too few arguments"
    exit 1
fi

# Compile the C++ file
g++ -L/opt/intel/mkl/lib/intel64 -I/opt/intel/mkl/include IntelBLAS.cpp -o IntelBLAS -fopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl

# Check if the compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "InelBLAS with $2 THREADS:"
# Run the executable 10 times
for i in $( eval echo {1..$1} ); do
    echo "Run #$i:"
    export MKL_NUM_THREADS="$2"; ./IntelBLAS "$3" "$4" "$5"
    echo "------------------------"
done

rm IntelBLAS
