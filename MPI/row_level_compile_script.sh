#!/bin/bash

# Check argument
if [ $# -lt 5 ]; then
  echo "Too few arguments"
    exit 1
fi

# Compile the C++ file
mpicxx row_level.cpp -o row_level

# Check if the compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "row_level with $2 THREADS:"
for i in $( eval echo {1..$1} ); do
    echo "Run #$i:"
    mpirun -np "$2" ./row_level "$3" "$4" "$5"
    echo "------------------------"
done

rm row_level