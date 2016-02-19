do for [t=1:100]{
file = sprintf('../Results/uGC%d.out',t)
tit = sprintf('t=%d',t)
set title tit
plot file u 2:3, "../Results/uGC_ref.out" u 1:2
pause(.1)
}
