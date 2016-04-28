#Plotting $F_z$ vs $n$ for the T=7 virus-indentation data
# $F_z$ is the force in z-direction no the Atomic Force 
# Microscopy indenter and $n$ is the step number

#Setting up the terminal
set terminal postscript enhanced eps color;
set output 'Fz.eps';
set size 1,1

# Setting up the axes
set title "Force vs Normalized Indentation"
set xlabel 'Normalized Indentation, z'
set ylabel 'Force, {/Times-Italic F}_z, (N)'            # italics!
set xtics 0,0.1,1
set mxtics 5
set ytics 0,0.1,0.4
set mytics 2
set grid xtics mxtics ytics mytics

#Now plotting
plot "Loading.fz" using 1:3 title "Force on AFM while loading" with linespoints pt 1 pi 0,\
     "Unloading.fz" using 1:3 title "Force on AFM while unloading" with linespoints pt 1 pi 0;


