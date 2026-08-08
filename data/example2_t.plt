set encoding iso_8859_1
set term postscript eps enhanced "Times-Roman" 20
set output "example2_t.eps"
set size 0.82
set logscale y
set key spacing 3
set xrange [0:0.2]
set yrange [1e-16:1]
set format y "10^{%L}"
set xlabel "{/Times-Roman=24 computation time [s]}"
set ylabel "{/Times-Roman=24 maximum error}"
plot "RZ_coll_ex2.dat" using 3:2 w lp title "SE-collocation (RZ)" lt 3 pt 8 ps 1.7, "SE_coll_ex2.dat" using 3:2 w lp title "SE-collocation (Stenger)" lt 3 pt 4 ps 1.7, "DE_coll_ex2.dat" using 3:2 w lp title "DE-collocation" lt 4 pt 6 ps 1.7, "SE_nyst_ex2.dat" using 3:2 w lp title "SE-Nystr{\366}m" lt 3 pt 2 ps 1.7, "DE_nyst_ex2.dat" using 3:2 w lp title "DE-Nystr{\366}m" lt 4 pt 1 ps 1.7
