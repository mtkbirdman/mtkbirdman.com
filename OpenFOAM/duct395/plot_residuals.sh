#!/bin/bash

#export DISPLAY=localhost:0.0
gnuplot -persist > /dev/null 2>&1 << EOF
set title "Residuals vs. Time"
set xlabel "Time / Iteration"
set ylabel "residuals [-]"
set logscale y
plot "plot_residuals.txt" using 1:2 title 'U_x' with lines,"plot_residuals.txt" using 1:3 title 'U_y' with lines,"plot_residuals.txt" using 1:4 title 'k' with lines,"plot_residuals.txt" using 1:5 title 'p' with lines
EOF