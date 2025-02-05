set terminal wxt
set xlabel "x"
set ylabel "u"
set title "Log plot of u vs N"
set grid

plot "error_J.dat" using 1:2 with lines title "error_J", \
     "error_GS.dat" using 1:2 with lines title "error_{GS}"

pause -1
