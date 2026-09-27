# Carte de B_r photosphérique, à partir de la sortie binaire de gong_export
# (float64, 360 x 180, longitude la plus rapide, sinus de latitude du sud au nord).
# Usage :
#   gnuplot -e "in='data/br_2232_129.bin'; out='figures/br_2232_129.png'; ttl='...'" plots/br_map.gp
NLON = 360
NLAT = 180
BMAX = 20   # G : saturation de l'échelle, pour voir les champs faibles (pôles)

set terminal pngcairo size 1400,720 font "Helvetica,13"
set output out
set title ttl
set xlabel "Carrington longitude (°)"
set ylabel "Latitude (°)"
set xrange [0:360]
set yrange [-1:1]
set xtics 0,60,360 out nomirror
set ytics out nomirror ("-90" -1, "-60" -sqrt(3)/2, "-30" -0.5, "0" 0, "30" 0.5, "60" sqrt(3)/2, "90" 1)
set border lc rgb "#6b6b68"
set cbrange [-BMAX:BMAX]
set cblabel "B_r (G)   blue: inward   red: outward"
set palette defined (-1 "#2a78d6", 0 "#f0efec", 1 "#e34948")
unset key

plot in binary array=(NLON,NLAT) format="%float64" \
     dx=360.0/NLON dy=2.0/NLAT origin=(180.0/NLON, -1 + 1.0/NLAT) with image
