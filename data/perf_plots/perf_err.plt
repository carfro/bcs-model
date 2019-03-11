# Plotting script for the performance benchmark of the Pfaffians
# Plot by loading in gnuplot then typeset in pdflatex;
# gnuplot load 'perf_err.plt' &&  pfdlatex perf_err && pdflatex perf_time

reset

# Epslatex
set terminal epslatex size 24cm,18cm color colortext standalone 

# Error/Precision plot
set output "perf_err.tex"
set title "Precision as {\\tiny $\\prod\\limits_{i}^{N/2} v_i^2:=\\text{Prod}(v_i)$}, {\\tiny $| \\frac{\\text{Prod}(v_i) - \\text{Pf}_x (W^T W)}{\\text{Prod}(v_i)} |$ }" font ",20"
#set title "Precision as {\\tiny $| \\frac{\\prod\\limits_{i}^{N/2} v_i^2 - \\text{Pf}_x (W^T W)}{\\prod\\limits_{i}^{N/2} v_i^2 } |$ }" font ",20"

# Color definitions
#set border linewidth 1.5
#set style line 1 lc rgb '#800000' lt 1 lw 2
#set style line 2 lc rgb '#ff0000' lt 1 lw 2
#set style line 3 lc rgb '#ff4500' lt 1 lw 2

set key outside
#set nokey

# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 1
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

#set mxtics 2
#set xrange [0.6:1.9]
#set xtics 0.1,0.2,1.9
set xlabel 'Size of matrix $n:W_{n \times n}$'
set ylabel 'Relative error' offset 0,0

set ytics out offset -1.25,0
set yrange [1e-18:1e-12]
set logscale y
set format x '{\small$%g$}'
set format y '{\small$10^{%T}$}'
set mytics 10

set key autotitle columnhead
plot for [i=3:5] 'performance_err.dat' u 1:i  
#plot for [i=3:5] 'performance_err.dat' u 1:i  title "Pf$_{".(i-2)."}$"
#plot for [i=3:5] 'performance_err.dat' u 1:( abs( abs(column(i))- $2 ) ) title "Pf$_{".(i-2)."}$"
#plot for [i=3:5] 'performance_err.dat' u 1:( abs(column(i))  )

set output

# Runtime plot --------------------------------------------------------------------------
set output "perf_time.tex"
set title "Runtime of Pfaffian routines" font ",20"

# Color definitions
#set border linewidth 1.5
#set style line 1 lc rgb '#800000' lt 1 lw 2
#set style line 2 lc rgb '#ff0000' lt 1 lw 2
#set style line 3 lc rgb '#ff4500' lt 1 lw 2

set key outside

# Axes
set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror out scale 0.75
# Grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

#set mxtics 2
#set mytics 10

#set xrange [0.6:1.9]
#set xtics 0.1,0.2,1.9
set xlabel 'Size of matrix $n:W_{n \times n}$'
set ylabel 'Runtime / s ' offset 1,-0.5
set logscale y
set ytics out offset -1.25,0
set yrange [1e-5:1e+5]
set format '{\small$%g$}'
set format y '{\small$10^{%T}$}'


plot for [i=3:5] 'performance_time.dat' u 1:( abs( abs(column(i))- $2 ) ) title "Pf$_{".(i-2)."}$"
#plot for [i=3:5] 'performance_err.dat' u 1:( abs(column(i))  )

set output















# ----------------------------------------------------------

#set terminal epslatex size 9cm,7cm color colortext standalone header \
#set output "perf_err.tex"
#set output "perf_err.ps"
#set size 1,1
#set title "Error, $|\prod_{i=1}^{N/2} v_i^2 - \text{Pf}(W^T W)$" font ",20"
#set key outside
#set format y '%g'
#set logscale y
#plot for [i=3:5] 'performance_err.dat' u 1:( abs( abs(column(i))- $2 ) ) title "bla"
#plot for [i=3:5] 'performance_err.dat' u 1:( abs(column(i))  )


#do for [i=2:5] {print '$'.i  }
#plot for [i=2:5] 'performance_err.dat' u 1:(x=sprintf("$%d",i) , abs(x))
#plot for [i=2:5] 'performance_err.dat' u 1:(abs(sprintf("$%d",i)))
#plot 'performance_err.dat' u 1:3
#replot 'performance_err.dat' u 1:(abs($3))
