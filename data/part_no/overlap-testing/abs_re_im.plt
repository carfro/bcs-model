reset
set term x11 1
plot 'overlap_imag.dat' u 2:3 w l title 'imag(pfaffian-pfa)', 'overlap_imag.dat' u 2:5 w l title 'imag(Pfaffian-ext)','overlap_imag.dat' u 2:4 w l title 'imag(Onishi)'
set term x11 2
plot 'overlap_real.dat' u 2:3 w l title 'real(pfaffian-pfa)', 'overlap_real.dat' u 2:5 w l title 'real(Pfaffian-ext)','overlap_real.dat' u 2:4 w l title 'real(Onishi)'
set term x11 3
plot 'overlap_proj_test.dat' u 2:3 w l title 'abs(pfaffian-pfa)', 'overlap_proj_test.dat' u 2:5 w l title 'abs(Pfaffian-ext)','overlap_proj_test.dat' u 2:4 w l title 'abs(Onishi)'

