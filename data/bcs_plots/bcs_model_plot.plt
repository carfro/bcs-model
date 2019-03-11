# Plotting script for the performance benchmark of the Pfaffians
# Plot by loading in gnuplot then typeset in pdflatex;
# gnuplot load 'perf_err.plt' &&  pfdlatex perf_err && pdflatex perf_time

reset

# Epslatex
#set terminal epslatex size 24cm,18cm color colortext standalone 

# Error/Precision plot
#set output "perf_err.tex"
#set title "Precision as {\\tiny $\\prod\\limits_{i}^{N/2} v_i^2:=\\text{Prod}(v_i)$}, {\\tiny $| \\frac{\\text{Prod}(v_i) - \\text{Pf}_x (W^T W)}{\\text{Prod}(v_i)} |$ }" font ",20"

# Color definitions
#set border linewidth 1.5
#set style line 1 lc rgb '#800000' lt 1 lw 2
#set style line 2 lc rgb '#ff0000' lt 1 lw 2
#set style line 3 lc rgb '#ff4500' lt 1 lw 2

#set key outside
#set nokey

# Axes
#set style line 11 lc rgb '#808080' lt 1
#set border 3 back ls 11
#set tics nomirror out scale 1
# Grid
#set style line 12 lc rgb'#808080' lt 0 lw 1
#set grid back ls 12

#set mxtics 2
set xrange [11:111]
#set xtics 0.1,0.2,1.9
#set xlabel 'Size of matrix $n:W_{n \times n}$'
#set ylabel 'Relative error' offset 0,0

#set ytics out offset -1.25,0
set ytics 0.1
set yrange [0:1]
#set format x '{\small$%g$}'
#set format y '{\small$10^{%T}$}'
#set mytics 10

set key autotitle columnhead
plot for [IDX=1:7] 'occupation.dat' i IDX u 1:2 w l title columnheader(1)  

#set output

