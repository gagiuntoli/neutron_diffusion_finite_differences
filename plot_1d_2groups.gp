# gnuplot -e "datafile='my_data.txt'" plot_3d.gp

set terminal pngcairo size 800,600
set output 'output_1d_2g.png'

# Titles and labels
set title "Comparison of val1 and val2"
set xlabel "X"
set ylabel "Values"

# Grid for better readability
set grid

# Legend positioning
set key top left

# Plot data
plot datafile using 1:2 with linespoints title "phi 1", \
     datafile using 1:3 with linespoints title "phi 2"