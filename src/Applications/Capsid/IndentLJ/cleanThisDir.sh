#!/bin/bash

# This script deletes the files created as a result of running a
# simulation
# -Amit

# The following two variables store the directory name and the list of
# files to be removed
DIRECTORY=.
REMOVEFILES=('abnormal*.vtk' 'contact*.vtk' 'lbfgsb*.vtk' \
'*-body*step*.vtk' '*.fz' '*.initial*.vtk' '*.relaxed*.vtk' \
'*~')

# We are passing all elements of $REMOVEFILES array and piping the
# output of ls to null device to make its execution silent
echo 'Following files will be removed:'
ls ${REMOVEFILES[@]} 2>/dev/null

echo 'Do you want to delete these files?[y/n]'
read deleteFlag1

if [ "$deleteFlag1" == "y" ]
then
    rm ${REMOVEFILES[@]} 2>/dev/null
    echo 'Directory is all cleaned up and ready for next simulation!'
else
    echo 'Clean up aborted. Call me again when you need it! See ya...'
fi