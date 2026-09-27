# Courbes en L (fichiers *_lcurve.txt de fit_map) : norme du résidu |A c - b| en
# abscisse, norme de la solution |c| en ordonnée ; lambda croît vers la droite.
# Usage :
#   gnuplot -e "full='data/cr2233_full_lcurve.txt'; cfull=K1; nopoles='data/cr2233_nopoles_lcurve.txt'; cnopoles=K2; out='figures/lcurve_cr2233.png'" plots/lcurve.gp
# K1, K2 : indice du coin, écrit en tête de chaque fichier (« # coin K »).
set terminal pngcairo size 1000,700 font "Helvetica,13"
set output out
set title "L-curves, GONG CR 2233, l_{max} = 30 (dot: corner = chosen {/Symbol l})"
set logscale xy
set autoscale xfix
set format x "%g"
set format y "10^{%L}"
set xlabel "Residual norm |Ac - b| (G)"
set ylabel "Solution norm |c| (G)"
set grid lc rgb "#e5e4e0"
set border lc rgb "#6b6b68"
set key top right

plot full u 2:3 w l lw 2 lc rgb "#2a78d6" t "full map (all latitudes)", \
     full u 2:3 every ::cfull::cfull w p pt 7 ps 2 lc rgb "#2a78d6" notitle, \
     nopoles u 2:3 w l lw 2 lc rgb "#eb6834" t "poles masked (latitudes above 60° excluded)", \
     nopoles u 2:3 every ::cnopoles::cnopoles w p pt 7 ps 2 lc rgb "#eb6834" notitle
