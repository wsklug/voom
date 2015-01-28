#Plotting $F_z$ vs $n$ for the T=7 virus-indentation data
# $F_z$ is the force in z-direction no the Atomic Force 
# Microscopy indenter and $n$ is the step number

#Setting up the terminal
set terminal postscript enhanced color;
set output 'Fz.ps';
set size 1,1

# Setting up the axes
set title "Force versus Step Number"
set xlabel 'Step, n'
set ylabel 'Force, {/Times-Italic F}_z, (N)'            # italics!
set xtics 0,100,1000
set mxtics 2
set ytics 0,0.1,0.4
set mytics 2
set grid xtics mxtics ytics mytics

#Now plotting
plot "T7input.fz" using 4:3 title 'Force on AFM' with linespoints pt 1 pi 0;

#Post-processing outside gnuplot
#latex Fz.tex
#dvips Fz.dvi
#ps2eps Fz.ps
