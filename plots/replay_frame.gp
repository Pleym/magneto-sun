# Une image de l'animation de la campagne (scripts/make_replay.sh) : B_r du modèle PFSS
# sur la surface source et sa ligne neutre, base de la nappe de courant héliosphérique.
# Variables : out, ss (carte *_ss.bin), neutral (*_neutral.txt), date, cr
NLON = 360
NLAT = 180

set terminal pngcairo size 1280,720 font "Sans,14"
set output out

set title sprintf("Solar wind source surface   %s   (Carrington rotation %s)", date, cr) font "Sans,18"
set tmargin at screen 0.90
set bmargin at screen 0.10
set lmargin at screen 0.09
set rmargin at screen 0.85
set xrange [0:360]
set yrange [-1:1]
set xtics 0,60,360 out nomirror
set ytics out nomirror ("-60°" -sqrt(3)/2, "-30°" -0.5, "0°" 0, "30°" 0.5, "60°" sqrt(3)/2)
set xlabel "Carrington longitude (°)"
set ylabel "Latitude"
set cbrange [-0.15:0.15]
set cbtics 0.05
set cblabel "Radial magnetic field (G)"
set palette defined (-1 "#2a78d6", 0 "#f0efec", 1 "#e34948")
unset key

plot ss binary array=(NLON,NLAT) format="%float64" dx=360.0/NLON dy=2.0/NLAT \
         origin=(180.0/NLON, -1 + 1.0/NLAT) with image notitle, \
     neutral u 1:2 w p pt 7 ps 0.55 lc rgb "black" notitle
