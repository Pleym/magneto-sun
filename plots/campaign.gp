# Produit L4 : accord de polarité par fenêtre sur toute une campagne, avec la position de
# la sonde (fichier de scripts/campaign_summary.sh). Colonnes : 1 fenêtre, 2 statut,
# 3 début, 6 accord de secteur, 7 référence triviale, 10-11 distance, 12-13 latitude.
# Usage : gnuplot -e "in='products/campaign_summary.txt'; out='...'" plots/campaign.gp
set terminal pngcairo size 1400,950 font "Helvetica,13"
set output out
set multiplot layout 3,1 title "Solar Orbiter vs PFSS (GONG): sector polarity agreement per 27-day window"
set xdata time
set timefmt "%Y-%m-%dT%H:%M"
set format x "%b %Y"
set lmargin 14
set rmargin 4
set grid lc rgb "#e5e4e0"
set border lc rgb "#6b6b68"
isOk = 'strcol(2) eq "ok"'

set ylabel "Agreement (%)"
set yrange [0:100]
set key bottom right
plot in u 3:(@isOk ? $6 : NaN) w lp lw 2 pt 7 lc rgb "#2a78d6" t "PFSS prediction", \
     in u 3:(@isOk ? $7 : NaN) w lp lw 2 pt 6 dt 2 lc rgb "#6b6b68" t "majority-polarity baseline"

set autoscale y
unset key
set ylabel "Latitude (°)"
plot in u 3:(@isOk ? ($12 + $13) / 2 : NaN):12:13 w yerrorbars pt 7 lc rgb "#1baf7a" notitle

set ylabel "Distance (AU)"
plot in u 3:(@isOk ? ($10 + $11) / 2 : NaN):10:11 w yerrorbars pt 7 lc rgb "#eb6834" notitle

unset multiplot
