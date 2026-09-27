# Polarité mesurée par Solar Orbiter et prédite par le PFSS (sortie de polarity).
# Colonnes : 1 temps, 6 v_R, 11 B le long de la spirale, 12 polarité horaire,
# 13 polarité de secteur, 15 polarité prédite.
# Usage : gnuplot -e "in='data/solo_series_full.txt'; out='...'; ttl='...'" plots/polarity_series.gp
OUTWARD = 0xe34948
INWARD = 0x2a78d6
HALF_HOUR = 1800

set terminal pngcairo size 1400,950 font "Helvetica,13"
set output out
set multiplot layout 3,1 title ttl
set xdata time
set timefmt "%s"
set format x "%d %b"
set autoscale xfix
set lmargin 16
set rmargin 4
set grid xtics lc rgb "#e5e4e0"
set border lc rgb "#6b6b68"

set ylabel "B along Parker spiral (nT)"
plot in u 1:($12 != 0 ? $11 : NaN):($11 > 0 ? OUTWARD : INWARD) w impulses lw 2 lc rgb variable notitle, \
     0 w l lc rgb "#6b6b68" notitle

set ylabel ""
set yrange [-0.6:1.6]
set ytics ("PFSS prediction" 0, "measured sector" 1)
set key outside top center horizontal
plot in u 1:($13 != 0 ? 1 : NaN):(HALF_HOUR):(0.35):($13 > 0 ? OUTWARD : INWARD) \
         w boxxyerror fs solid noborder lc rgb variable notitle, \
     in u 1:($15 != 0 ? 0 : NaN):(HALF_HOUR):(0.35):($15 > 0 ? OUTWARD : INWARD) \
         w boxxyerror fs solid noborder lc rgb variable notitle, \
     NaN w boxes fs solid lc rgb OUTWARD t "outward (B_r > 0)", \
     NaN w boxes fs solid lc rgb INWARD t "inward (B_r < 0)"

set autoscale y
set ytics auto
unset key
set ylabel "Solar wind v_R (km/s)"
plot in u 1:6 w l lw 2 lc rgb "#1a1a19" notitle

unset multiplot
