#!/bin/bash

# Check argument
if [ $# -lt 5 ]; then
  echo "Too few arguments"
    exit 1
fi

# Compile the C++ file
g++ -o collapse_three  collapse_three.cpp -fopenmp

# Check if the compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "Collapse three with $2 THREADS:"
# Run the executable 10 times
for i in $( eval echo {1..$1} ); do
    echo "Run #$i:"
    export OMP_NUM_THREADS="$2"; ./collapse_three "$3" "$4" "$5"
    echo "------------------------"
done

rm collapse_three