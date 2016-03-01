Re=100

set view map
set pm3d
unset surface

set xrange[0:1]
set yrange[0:1]
set size square

set contour 
set cntrparam levels discrete -1e-10, -1e-7, -1e-5, -1e-4, -.01, -.03, -.05,\
-.07, -.09, -.1, -.11, -.115, -.1175, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4,\
2.5e-4, 5e-4, 1e-3, 1.5e-3, 3e-3

splot sprintf("../Results/psi_N129_Re%d.out", Re) u 1:2:3 w lines lc rgb "black",\
    sprintf("../Results/vec_N129_Re%d.out", Re) every 50:50 \
    u 1:2:(0):($3*0.10):($4*0.10):(0) w vectors lc rgb "white"
pause(-1)
