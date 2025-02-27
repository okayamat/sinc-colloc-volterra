set encoding iso_8859_1
set term postscript eps enhanced "Times-Roman" 20
set output "example1_N.eps"
set size 0.82
set logscale y
set key spacing 3
set xrange [0:200]
set yrange [1e-16:10]
set format y "10^{%L}"
set xlabel "{/Times-Italic=24 N}"
set ylabel "{/Times-Roman=24 maximum error}"
plot "RZ_coll_ex1.dat" using 1:2 w lp title "SE-Sinc-collocation (Rashidinia-Zarebnia)" lt 3 pt 8 ps 1.7, "SE_coll_ex1.dat" using 1:2 w lp title "SE-Sinc-collocation (Stenger)" lt 3 pt 4 ps 1.7, "DE_coll_ex1.dat" using 1:2 w lp title "DE-Sinc-collocation" lt 4 pt 6 ps 1.7, "SE_nyst_ex1.dat" using 1:2 w lp title "SE-Sinc-Nystr{\366}m" lt 3 pt 2 ps 1.7, "DE_nyst_ex1.dat" using 1:2 w lp title "DE-Sinc-Nystr{\366}m" lt 4 pt 1 ps 1.7
