#!/bin/bash

if [ $# -lt 5 ]; then
  echo "Too few arguments"
    exit 1
fi

# Loop through all .sh files in the current directory
for script in *.sh; do
  # Skip the script itself
  if [[ "$script" == "$0" ]]; then
    continue
  fi

  # Execute the script and pass along any arguments
  echo "Executing $script..."
  bash "$script" "$1" "$2" "$3" "$4" "$5"
done