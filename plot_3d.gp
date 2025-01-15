# gnuplot -e "datafile='my_data.txt'" plot_3d.gp

set terminal pngcairo enhanced size 800,600
set output 'output.png'

set xlabel "X-axis"
set ylabel "Y-axis"
set zlabel "Value"

set title "3D Plot of x, y, val"
set grid

splot datafile using 1:2:3 with points pointtype 7 pointsize 1 palette notitle
