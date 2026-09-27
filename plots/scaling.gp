# Passage à l'échelle forte : accélération T(1) / T(n) selon le nombre de threads
# (CSV de scripts/bench_scaling.sh) : anneaux sur HMI (lmax 720, OpenMP) et dense sur
# GONG (lmax 60, BLAS + OpenMP). En option, une deuxième série d'anneaux (in2, légende
# lab2), par exemple la même machine avec une autre BLAS.
# Colonnes : 1 carte, 2 threads BLAS, 3 solveur, 4 lmax, 6 threads OpenMP, 11 total.
# Usage : gnuplot -e "in='data/bench_mac_m3.csv'; lab='...'; out='...'; ttl='...'" plots/scaling.gp
#   (ajouter in2='...'; lab2='...' pour la deuxième série)
set datafile separator ","
if (!exists("in2")) in2 = ''
if (!exists("lab")) lab = ''
if (!exists("lab2")) lab2 = ''
# Une série = un cas ; smooth unique moyenne les répétitions
rings(file) = sprintf("< awk -F, '$3 == \"rings\" && $4 == 720' %s", file)
dense(file) = sprintf("< awk -F, '$3 == \"dense\" && $4 == 60 && $2 == $6' %s", file)
stats rings(in) u ($6 == 1 ? $11 : NaN) nooutput
t1Rings = STATS_mean
stats dense(in) u ($6 == 1 ? $11 : NaN) nooutput
t1Dense = STATS_mean
stats rings(in) u 6 nooutput
maxThreads = STATS_max
t1Rings2 = 1
if (in2 ne '') {
    stats rings(in2) u ($6 == 1 ? $11 : NaN) nooutput
    t1Rings2 = STATS_mean
}

set terminal pngcairo size 1000,700 font "Helvetica,13"
set output out
set title ttl
set xlabel "Threads"
set ylabel "Speedup T(1) / T(n)"
set xrange [0.8:maxThreads * 1.1]
set yrange [0:maxThreads * 1.1]
set grid lc rgb "#e5e4e0"
set border lc rgb "#6b6b68"
set key top left

plot x w l dt 2 lc rgb "#6b6b68" t "ideal", \
     rings(in) u 6:(t1Rings / $11) smooth unique w lp lw 2 pt 7 lc rgb "#1baf7a" \
         t sprintf("rings, HMI l_{max} = 720, %s", lab), \
     dense(in) u 6:(t1Dense / $11) smooth unique w lp lw 2 pt 7 lc rgb "#2a78d6" \
         t sprintf("dense, GONG l_{max} = 60, %s", lab), \
     (in2 ne '' ? rings(in2) : rings(in)) u 6:(in2 ne '' ? t1Rings2 / $11 : NaN) smooth unique \
         w lp lw 2 pt 5 lc rgb "#eb6834" t (in2 ne '' ? sprintf("rings, HMI l_{max} = 720, %s", lab2) : "")
