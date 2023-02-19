reset

set terminal jpeg enhanced size 800,600 font "Arial,13"
set output 'conveccionLTS.jpeg'

set title '' 
set xlabel 'x' 
set ylabel 'v' 

#set logscale x

set xtics nomirror
set ytics nomirror
#set format y "%.1f"
#set format x "%.1f"
#set grid

set key top left

set xrange [0:80]
set yrange [0:6]

#f(x) = x <= 25 ? 0. : (x <= 55 ? 5. : 0.)
f(x) = 5.* exp(-(x-40)**2/2./5.**2)

plot f(x-10) t 'Analytic' lc 'black', 'OutputFiles/conveccionLinealLocalTimeStep1.dat' u 1:2 t 'LTS1' pt 7 ps 1 lc 'red', 'OutputFiles/conveccionLinealUpwind1.dat' u 1:2 t 'Upwind' pt 7 ps 1 lc 'web-blue'