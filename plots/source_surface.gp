# B_r du modèle PFSS sur la surface source (sortie *_ss.bin de pfss_map), sa ligne
# neutre, et en option la trace des pieds balistiques de la sonde (Phase 4),
# colorée par la polarité de secteur mesurée (fichier *_series.txt de polarity).
# Usage :
#   gnuplot -e "in='data/cr2233_full_ss.bin'; neutral='data/cr2233_full_neutral.txt'; out='...'; ttl='...'; bmax=0.15" plots/source_surface.gp
#   (ajouter track='data/solo_series.txt' pour la trace de la sonde)
NLON = 360
NLAT = 180
if (!exists("track")) track = ''
OUTWARD = 0xe34948
INWARD = 0x2a78d6

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
set cbrange [-bmax:bmax]
set cblabel "B_r at 2.5 R_{sun} (G)"
set palette defined (-1 "#2a78d6", 0 "#f0efec", 1 "#e34948")
set key bottom right opaque box

map = sprintf("%s", in)
if (track eq '') {
    plot map binary array=(NLON,NLAT) format="%float64" dx=360.0/NLON dy=2.0/NLAT \
             origin=(180.0/NLON, -1 + 1.0/NLAT) with image notitle, \
         neutral u 1:2 w p pt 7 ps 0.35 lc rgb "#1a1a19" t "neutral line (base of the HCS)"
} else {
    # Colonnes de la série : 9 = longitude du pied (°), 10 = latitude (°), 13 = polarité de secteur
    plot map binary array=(NLON,NLAT) format="%float64" dx=360.0/NLON dy=2.0/NLAT \
             origin=(180.0/NLON, -1 + 1.0/NLAT) with image notitle, \
         neutral u 1:2 w p pt 7 ps 0.35 lc rgb "#1a1a19" t "neutral line (base of the HCS)", \
         track u 9:($13 != 0 ? sin($10*pi/180) : NaN) w p pt 7 ps 1.1 lc rgb "#1a1a19" notitle, \
         track u 9:($13 != 0 ? sin($10*pi/180) : NaN):($13 > 0 ? OUTWARD : INWARD) w p pt 7 ps 0.7 lc rgb variable \
             t "Solar Orbiter footpoints, colored by measured polarity"
}
