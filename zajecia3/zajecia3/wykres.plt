set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'interpolacja.png'

set title "Interpolacja Lagrange'a – Zjawisko Rungego"
set xlabel "x"
set ylabel "y"
set grid

plot "dane.txt" using 1:2 with lines lw 2 title "f(x) = 1 / (1 + x^3)", \
     "rown.txt" using 1:2 with lines lw 2 lc rgb "red" title "Interpolacja – węzły równoodległe", \
     "czeb.txt" using 1:2 with lines lw 2 lc rgb "blue" title "Interpolacja – węzły Czebyszewa"
