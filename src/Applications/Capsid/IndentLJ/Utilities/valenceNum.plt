#Setting up the terminal
set terminal postscript enhanced eps color;
set output 'Defects.eps';
set size 1,1

# Setting up the axes
set title "Number of Defects vs Step Number"
set xlabel 'Step, n'
set ylabel 'Number of Defects, (N)'
#set xrange [0:2400]
set yrange [0:25]
set xtics 0,200,2400
set mxtics 2
#set ytics 0,50,300
set ytics 0,5,25
set mytics 5
set grid xtics mxtics ytics mytics

#By default gnuplot expects tab or space separated values in the file
#We need to specify that our file has comma-separated values
set datafile separator ","

#Now plotting
plot "valenceCount.csv" using 2 title 'Pentamers' with lines lw 2,\
	"valenceCount.csv" using 4 title 'Heptamers' with lines lw 2,\
	"valenceCount.csv" using 5 title 'Octamers' with lines lw 2,\
	"valenceCount.csv" using 6 title 'Net Charge' with lines lw 2;


