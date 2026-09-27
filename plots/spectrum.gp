# Spectre de B_r par degré l (fichier *_spectrum.txt de pfss_map) : à la surface et
# sur la surface source. Montre le filtrage passe-bas du modèle PFSS.
# Usage : gnuplot -e "in='data/cr2233_full_spectrum.txt'; out='figures/spectrum_cr2233.png'" plots/spectrum.gp
set terminal pngcairo size 1000,650 font "Helvetica,13"
set output out
set title "Mean-square B_r per degree l, GONG CR 2233 (full-map fit)"
set logscale y
set format y "10^{%L}"
set xlabel "Spherical harmonic degree l"
set ylabel "<B_r^2>_l (G^2)"
set xrange [0.5:*]
set grid lc rgb "#e5e4e0"
set border lc rgb "#6b6b68"
set key top right

plot in u 1:2 w lp lw 2 pt 7 ps 0.8 lc rgb "#2a78d6" t "photosphere (r = R_{sun})", \
     in u 1:3 w lp lw 2 pt 7 ps 0.8 lc rgb "#eb6834" t "source surface (r = 2.5 R_{sun})"
