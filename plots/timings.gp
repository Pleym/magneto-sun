# Temps de l'ajustement selon lmax, toute la machine (CSV de scripts/bench_scaling.sh) :
# solveur dense sur GONG, solveur par anneaux sur GONG et sur HMI.
# Colonnes : 1 carte, 2 threads BLAS, 3 solveur, 4 lmax, 5 pixels, 6 threads OpenMP, 11 total.
# Usage : gnuplot -e "in='data/bench_mac_m3.csv'; maxt=8; out='...'; ttl='...'" plots/timings.gp
set datafile separator ","
set terminal pngcairo size 1000,700 font "Helvetica,13"
set output out
set title ttl
set logscale xy
set format y "10^{%L}"
set xlabel "l_{max}"
set ylabel "Fit wall time (s)"
set xrange [8:1000]
set grid lc rgb "#e5e4e0"
set border lc rgb "#6b6b68"
set key top left

# Une série = une carte, un solveur, maxt threads ; smooth unique moyenne les répétitions
sel(map, solver) = sprintf("< awk -F, '$1 ~ /^%s/ && $3 == \"%s\" && $6 == %d' %s", map, solver, maxt, in)
plot sel("mrzqs", "dense") u 4:11 smooth unique w lp lw 2 pt 7 lc rgb "#2a78d6" t "dense, GONG 360 x 180", \
     sel("mrzqs", "rings") u 4:11 smooth unique w lp lw 2 pt 7 lc rgb "#eb6834" t "rings, GONG 360 x 180", \
     sel("hmi", "rings") u 4:11 smooth unique w lp lw 2 pt 7 lc rgb "#1baf7a" t "rings, HMI 3600 x 1440"
