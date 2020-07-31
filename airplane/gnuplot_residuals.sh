#!/bin/bash

export DISPLAY=localhost:0.0
gnuplot -persist > /dev/null 2>&1 << EOF
set terminal x11
set title "Residuals vs. Time"
set xlabel "Time / Iteration"
set ylabel "residuals [-]"
set logscale y
plot "gnuplot_residuals.txt" using 1:2 title 'U_x' with lines,"gnuplot_residuals.txt" using 1:3 title 'U_y' with lines,"gnuplot_residuals.txt" using 1:4 title 'U_z' with lines,"gnuplot_residuals.txt" using 1:5 title 'p' with lines
EOF