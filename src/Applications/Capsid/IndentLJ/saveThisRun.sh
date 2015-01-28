#!/bin/bash

# This script saves the files created as a result of running a
# simulation into ./GoodRuns
# -Amit

# The following two variables store the directory name and the list of
# files to be removed
DIRECTORY=.
SAVEFILES=('abnormal*.vtk' 'contact*.vtk' 'lbfgsb*.vtk' \
'*-body*step*.vtk' '*.fz' '*.initial*.vtk' '*.relaxed*.vtk' \
'indent.log')

# We are passing all elements of $SAVEFILES array and piping the
# output of ls to null device to make its execution silent
#echo 'Following files will be saved:'
#ls ${SAVEFILES[@]} 2>/dev/null

echo 'What folder to save these files in ./GoodRuns/?'
read FOLDERNAME

NEWFOLDER="./GoodRuns/$FOLDERNAME"

if [ ! -d "$folderName" ]
then
    mkdir -p "$NEWFOLDER"
    mv ${SAVEFILES[@]} "$NEWFOLDER" 2>/dev/null
    cd "$NEWFOLDER"
    tar -zcvf "$FOLDERNAME.tar.gz" *
    mv "$FOLDERNAME.tar.gz" ../
    echo 'Simulation archived successfully!'
else
    echo 'Something went wrong... try debugging or do it manually :-P'
fi