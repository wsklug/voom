#!/bin/bash

# This script saves the files created as a result of running a
# simulation into $SAVETO
# -Amit

# The following two variables store the directory name and the list of
# files to be removed

SAVEFILES=('abnormal*.vtk' 'contact*.vtk' 'lbfgsb*.vtk' \
'*-body*step*.vtk' '*.fz' '*.initial*.vtk' '*.relaxed*.vtk' \
'*.output.*' '*.joblog.*')

# We are passing all elements of $SAVEFILES array and piping the
# output of ls to null device to make its execution silent
#echo 'Following files will be saved:'
#ls ${SAVEFILES[@]} 2>/dev/null

# Current directory is the default source directory
defSrcDir=$PWD
srcDir=""${1:-$defSrcDir}
echo "Source directory: $srcDir"

defTarDir=~/Research/Results/VirusIndentation/${srcDir##*/}
target=${2:-$defTarDir}
echo "Target directory: ${target}"

mkdir -p $target

for i in ${SAVEFILES[@]}
do
    mv ${srcDir}/${i} "$target" 2>/dev/null         
done

cd "$target"
tar -zcvf "${target##*/}.tar.gz" *
mv "${target##*/}.tar.gz" ../
echo 'Simulation archived successfully!'