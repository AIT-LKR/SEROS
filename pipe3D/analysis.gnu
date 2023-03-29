set multiplot layout 6,2
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"

## iter
#rhobar: inlet/global
plot 'measurementsFlow0' using 2:4 with lines title columnheader
plot 'measurementsFlow0' using 2:5 with lines title columnheader

#kineticEnergy
plot 'measurementsFlow0' using 2:6 with lines title columnheader

#shearStress: min/max/diff
#plot 'measurementsFlow0' using 2:7 with lines title columnheader
#plot 'measurementsFlow0' using 2:8 with lines title columnheader
plot 'measurementsFlow0' using 2:9 with lines title columnheader

#hist: mean/sigma/skew/kurt
#plot 'measurementsFlow0' using 2:10 with lines title columnheader
#plot 'measurementsFlow0' using 2:11 with lines title columnheader
#plot 'measurementsFlow0' using 2:12 with lines title columnheader
#plot 'measurementsFlow0' using 2:13 with lines title columnheader

## seros
#rhobarInletAvg(Smooth)/Grad(Smooth)
plot 'measurementsSeros0' using 3:4 with lines title columnheader, 'measurementsSeros0' using 3:5 with lines title columnheader
plot 'measurementsSeros0' using 3:6 with lines title columnheader, 'measurementsSeros0' using 3:7 with lines title columnheader
#rhobarGlobalAvg(Smooth)/Grad(Smooth)
plot 'measurementsSeros0' using 3:8 with lines title columnheader, 'measurementsSeros0' using 3:9 with lines title columnheader
plot 'measurementsSeros0' using 3:10 with lines title columnheader, 'measurementsSeros0' using 3:11 with lines title columnheader

#kineticEnergy(Smooth)/Grad
plot 'measurementsSeros0' using 3:12 with lines title columnheader, 'measurementsSeros0' using 3:13 with lines title columnheader
plot 'measurementsSeros0' using 3:14 with lines title columnheader

#diameter
plot 'measurementsSeros0' using 3:15 with lines title columnheader
#fluid cells
plot 'measurementsSeros0' using 3:16 with lines title columnheader

pause 5
reread
