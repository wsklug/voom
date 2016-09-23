#Setting up the terminal
set terminal postscript enhanced color;
set output 'Defects.ps';
#set size 1,1

#Set two rows of sub-plots with under one column
set multiplot layout 2,1

# Setting up the axes
set title "Defects versus Step Number"
set xlabel 'Step, n'
set ylabel 'Number of each type, (N)'
#set xtics 0,100,1000
#set mxtics 2
#set ytics 0,10,220
#set mytics 2
#set grid xtics mxtics ytics mytics

#By default gnuplot expects tab or space separated values in the file
#We need to specify that our file has comma-separated values
set datafile separator ","

#Now plotting
plot "valenceCount.txt" using 1:2 title 'Pentamers' with lines lw 2,\
	"valenceCount.txt" using 1:3 title 'Hexamers' with lines lw 2,\
	"valenceCount.txt" using 1:4 title 'Heptamers' with lines lw 2,\
	"valenceCount.txt" using 1:5 title 'Octamers' with lines lw 2,\
	"valenceCount.txt" using 1:6 title 'Net Charge' with lines lw 2;

#Post-processing outside gnuplot
#latex Fz.tex
#dvips Fz.dvi
#ps2eps Fz.ps
