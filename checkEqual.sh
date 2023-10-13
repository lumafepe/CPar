#!/bin/bash

file1="cp_output.txt"
file2="output/cp_output.txt"

# Calculate MD5 checksums for the files
md5sum_file1=$(md5sum "$file1" | awk '{print $1}')
md5sum_file2=$(md5sum "$file2" | awk '{print $1}')

# Compare the MD5 checksums
if [ "$md5sum_file1" == "$md5sum_file2" ]; then
  echo "OH DEAR, OH DEAR, GORGEOUS."
else
  echo "YOU FUCKING DONKEY."
fi