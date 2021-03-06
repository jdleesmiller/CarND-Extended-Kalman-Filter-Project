#
# Run with:
# gnuplot -e "filename='data/output-1'" plot.gnuplot
# gnuplot -e "filename='data/output-2'" plot.gnuplot
#

print filename
datafile = filename.'.txt'
figurefile = filename.'.png'

set terminal png size 1024, 1024
set output figurefile

set multiplot
set size 1.0, 1.0
set origin 0.0, 0.0
set key outside

set title 'Estimated vs Ground Truth X and Y Positions'
set xlabel 'X Position'
set ylabel 'Y Position'
set origin 0.0, 0.5
set size 1.0, 0.5

plot datafile using 1:2 title 'est xy', datafile using 7:8 title 'gt xy'

set title 'Estimated vs Ground Truth X and Y Velocities'
set xlabel 'Time Index'
set ylabel 'X or Y Velocity'
set origin 0.0, 0.25
set size 1.0, 0.25

plot datafile using 3 with lines title 'est vx', \
  datafile using 9 with lines title 'gt vx', \
  datafile using 4 with lines title 'est vy', \
  datafile using 10 with lines title 'gt vy'

set title 'Estimated vs Ground Truth Speed'
set xlabel 'Time Index'
set ylabel 'Speed'
set origin 0.0, 0.0
set size 1.0, 0.25

plot datafile using (sqrt($3*$3+$4*$4)) with lines title 'est v ', \
  datafile using (sqrt($9*$9+$10*$10)) with lines title 'gt v'
