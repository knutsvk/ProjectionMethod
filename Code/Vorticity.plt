set view map
set pm3d
unset surface

set xrange[0:1]
set yrange[0:1]
set size square

set contour 
set cntrparam levels discrete -3,-2,-1,-.5,0,.5,1,2,3,4,5
splot "../Results/omega_N129_Re3200.out" u 1:2:3 with lines lc rgb "black"
pause(-1)
