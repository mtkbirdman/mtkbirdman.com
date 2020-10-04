#!/bin/bash

export DISPLAY=localhost:0.0
gnuplot -persist > /dev/null 2>&1 << EOF
set terminal x11
set title "Coeffs vs. Time"
set xlabel "Time / Iteration"
set ylabel "Coeffs [-]"
set yrange [-0.5:1.5]
plot "gnuplot_coeffs.txt" using 1:2 title 'Cl' with lines,"gnuplot_coeffs.txt" using 1:3 title 'Cd' with lines,"gnuplot_coeffs.txt" using 1:4 title 'Cm' with lines
EOF