set dgrid3d 30,30
set hidden3d

set view map
set pm3d
unset surface

splot "Results/test.out" u 1:2:3 with lines
pause(10)
