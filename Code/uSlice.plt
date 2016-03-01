#   n   Re
#   2   100
#   3   400
#   4   1000
#   5   3200

Re=3200
n=5

plot sprintf("../Results/uGC_N17_Re%d.out",Re) u 1:2,\
    sprintf("../Results/uGC_N33_Re%d.out",Re) u 1:2,\
    sprintf("../Results/uGC_N65_Re%d.out",Re) u 1:2,\
    sprintf("../Results/uGC_N129_Re%d.out",Re) u 1:2,\
    "../Results/uGC.ref" u 1:n,\
    sprintf("../Results/vGC_N17_Re%d.out",Re) u 1:2,\
    sprintf("../Results/vGC_N33_Re%d.out",Re) u 1:2,\
    sprintf("../Results/vGC_N65_Re%d.out",Re) u 1:2,\
    sprintf("../Results/vGC_N129_Re%d.out",Re) u 1:2,\
    "../Results/vGC.ref" u 1:n
pause(-1)
