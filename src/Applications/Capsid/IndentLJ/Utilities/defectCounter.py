#!/usr/bin/env python
#############################################################################
#This script counts the pentamers, hexamers and heptamers at every step of a#
#simulation. It takes as input the .vtk file name format.                   #
#                                                                           #
#Author: Amit Singh                                                         #
#Date: 01/13/2014                                                           #
#                                                                           #
#############################################################################

import re

#Generally, the VTK output files of my simulations have names like
#'T7input-body0-step0001.vtk' to 'T7input-body0-step2002.vtk'. We need
#to read each of these files. The valency data is written in the file
#in the following format:
#
#FIELD FieldData 1
#Valence 1 16 unsigned_int
#5 6 6 6 5 6
#6 6 6 6 6 6
#5 6 6 4
#
#There are 16 valence values as can be seen from the 3 word in the
#line starting with 'Valence'. Our logic is to 1. read the file till
#we reach the line that starts with 'Valence' 2. Extract the total
#number of valence values to read and then 3. from the next line
#classify the values as penta, hexa or heptamers till we have as many
#values (16 in this example) as expected

#The 'base word' of the file name is 'T7input-body0-step' in our example
fileBaseName = "T7input-body0-step"

#The step numbering starts from 0000 so actually there are 2003 files
#but we will enter 2002 below
steps = 2002

#We need to pad the step number with enough zeros e.g 0001, 0103 to
#get the required file name. We prepare the format specifier here to
#use it later in the for loop
digits = len(str(steps))
formatSpec = '{0:0' + str(digits) + 'd}'

#We prepare the output file to be written containing the tab-separated
#valence data
output = open("valenceCount.csv","w")
header = "#Step,Pentamer,Hexamer,Heptamer,Octamers,NetCharge"
print >> output, header

for step in range(0, steps+1):
    currFile = fileBaseName + formatSpec.format(step) + ".vtk"

    with open(currFile,"r") as file:
#Reset the flag that checks when to start reading valence data
        reached = False

#Reset all counters
        totalCount = 0
        pentamers = 0
        hexamers = 0
        heptamers = 0
        octamers = 0
        charge = 0

        for line in file:

#Skip to next line as long as we don't reach line starting with 'Valence'
            if re.match("Valence.+",line):

#Extract the third word from the 'Valence' line to know how many
#values to look for
                words = line.split(" ")
                totalCount = int(words[2])

#Set the 'reached' flag on to start reading numbers from next line
                reached = True
                continue

            if reached and totalCount > 0:
                nums = line.strip().split(" ")
                for num in nums:                    
                    if int(num) == 5:
                        pentamers += 1                        
                    elif int(num) == 6:
                        hexamers += 1
                    elif int(num) == 7:
                        heptamers += 1
                    elif int(num) == 8:
                        octamers += 1
                    else:
                        print currFile
                        print "Absurd data found:" + num + "?\n";
                    totalCount -= 1
                    charge = charge + (6-int(num))
        
    printString = formatSpec.format(step) +","+ str(pentamers) + "," + \
        str(hexamers) + "," + str(heptamers) +","+str(octamers) +","+ \
        str(charge)
    print >> output, printString

output.close()
