#!/bin/bash

# Check argument
if [ $# -lt 4 ]; then
  echo "Too few arguments"
    exit 1
fi

# Compile the C++ file
g++ -o sequential_moltiplication sequential_moltiplication.cpp -fopenmp

# Check if the compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

echo "sequential_transpose:"
# Run the executable 10 times
for i in $( eval echo {1..$1} ); do
    echo "Run #$i:"
    ./sequential_moltiplication "$2" "$3" "$4"
    echo "------------------------"
done

rm sequential_moltiplication