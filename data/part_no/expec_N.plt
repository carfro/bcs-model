reset
set term x11 1
plot 'Nexp.dat' u 1:2 w l title 'Part nbr - Pfaffian', 'Nexp.dat' u 1:3 w l title 'Part nbr - Onishi', 'Nexp.dat' u 1:4 w l title 'Part nbr - Pfa_ext'
set term x11 2
plot 'rot_ol_sum_100l_24p_4d.dat' u 1:2 w l title 'Rot overlap - Pfaffian', 'rot_ol_sum_100l_24p_4d.dat' u 1:3 w l title 'Rot overlap - Onishi', 'rot_ol_sum_100l_24p_4d.dat' u 1:4 w l title 'Rot overlap - Pfa_{ext}'

