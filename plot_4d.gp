# gnuplot -e "datafile='my_data.txt'" plot_4d.gp

set terminal pngcairo size 800,600 enhanced
set output 'output_4d.png'

set palette rgbformulae 33,13,10  # Define the color palette (can be customized)
set cblabel "Value"               # Label for the color bar
set colorbox                      # Enable the color bar
set xlabel "X-axis"
set ylabel "Y-axis"
set zlabel "Z-axis"

set view 60,30                    # Set 3D view angles
set pointsize 1.5                 # Set fixed point size

splot datafile using 1:2:3:4 with points pointtype 7 palette notitle