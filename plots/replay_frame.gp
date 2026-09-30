# Une image du rejeu d'exploitation (scripts/make_replay.sh) : bandeau de la chaîne
# L1 -> L4 avec l'état de chaque étape, carte de la surface source de la fenêtre en
# cours, et graphique de campagne qui se complète fenêtre après fenêtre.
# Variables : out, header, s1..s4 (done | active | wait | warn), d1..d4 (détails),
#   ss, neutral, track, maptitle (carte ; ss = '' : pas encore de carte),
#   summary (synthèse L4), nplot (fenêtres affichées), ncur (fenêtre mise en valeur, 0 : aucune)
NLON = 360
NLAT = 180
OUTWARD = 0xe34948
INWARD = 0x2a78d6
BLUE = "#2a78d6"
GREY = "#6b6b68"
INK = "#1a1a19"
stateColor(s) = s eq "done" ? "#0ca30c" : s eq "active" ? "#2a78d6" : s eq "warn" ? "#fab219" : "#b5b3ad"
stateText(s) = s eq "done" ? "done" : s eq "active" ? "running" : s eq "warn" ? "insufficient data" : "waiting"

set terminal pngcairo size 1280,720 font "Helvetica,15" background rgb "#fcfcfb"
set output out
set multiplot

# --- Bandeau : titre et chaîne L1 -> L4, sur un tracé vide plein cadre -----------------
set origin 0,0
set size 1,1
set margins 0,0,0,0
unset border
unset tics
unset key
set xrange [0:1]
set yrange [0:1]
set label 1 "magneto-sun · ground-segment-style processing chain" \
    at screen 0.025,0.955 font "Helvetica-Bold,19" tc rgb INK
set label 2 header at screen 0.975,0.955 right font "Helvetica,14" tc rgb GREY

array TITLE[4] = ["L1 · raw inputs", "L2 · standardised", "L3 · science products", "L4 · campaign"]
array SUB[4] = ["GONG map · MAG · SWA-PAS · orbit", "hourly in situ series", \
                "SH fit, PFSS, polarity score", "agreement per window"]
array STATE[4] = [s1, s2, s3, s4]
array DETAIL[4] = [d1, d2, d3, d4]
x0 = 0.025
w = 0.215
gap = 0.035
ybot = 0.75
ytop = 0.89
do for [i = 1:4] {
    xa = x0 + (i - 1) * (w + gap)
    xb = xa + w
    isActive = STATE[i] eq "active"
    textColor = isActive ? "#ffffff" : INK
    set object i rect from screen xa,ybot to screen xb,ytop fc rgb stateColor(STATE[i]) \
        fs solid (isActive ? 0.9 : 0.15) border lc rgb stateColor(STATE[i]) lw 2
    set label 10 + i TITLE[i] at screen xa + 0.01, ytop - 0.035 font "Helvetica-Bold,15" tc rgb textColor
    set label 20 + i SUB[i] at screen xa + 0.01, ytop - 0.072 font "Helvetica,11" tc rgb textColor
    set label 30 + i stateText(STATE[i]) at screen xa + 0.01, ybot + 0.022 font "Helvetica-Bold,12" tc rgb textColor
    set label 40 + i DETAIL[i] at screen xa + 0.004, ybot - 0.035 font "Helvetica,12" tc rgb GREY
    if (i < 4) {
        set arrow i from screen xb + 0.004, (ybot + ytop) / 2 to screen xb + gap - 0.004, (ybot + ytop) / 2 \
            lw 2 lc rgb GREY size screen 0.008,25 filled
    }
}
plot NaN notitle
unset object
unset label
unset arrow

# --- Carte L3 : surface source de la fenêtre en cours ---------------------------------
set border lc rgb GREY
set origin 0.0,0.02
set size 0.53,0.66
set lmargin 8
set rmargin 11
set tmargin 2
set bmargin 3
set xrange [0:360]
set yrange [-1:1]
set xtics 0,90,360 out nomirror font ",11"
set ytics out nomirror font ",11" ("-60°" -sqrt(3)/2, "-30°" -0.5, "0°" 0, "30°" 0.5, "60°" sqrt(3)/2)
set xlabel "Carrington longitude (°)" font ",12"
set title maptitle font ",13" tc rgb INK
set cbrange [-0.15:0.15]
set cbtics 0.1 font ",10"
set cblabel "B_r at 2.5 R_{sun} (G)" font ",11"
set palette defined (-1 "#2a78d6", 0 "#f0efec", 1 "#e34948")
if (ss eq '') {
    unset colorbox
    set label 50 "waiting for L3 products…" at graph 0.5,0.5 center tc rgb GREY
    plot NaN notitle
    unset label 50
} else {
    set colorbox
    plot ss binary array=(NLON,NLAT) format="%float64" dx=360.0/NLON dy=2.0/NLAT \
             origin=(180.0/NLON, -1 + 1.0/NLAT) with image notitle, \
         neutral u 1:2 w p pt 7 ps 0.25 lc rgb INK notitle, \
         track u 9:($13 != 0 ? sin($10 * pi / 180) : NaN) w p pt 7 ps 0.9 lc rgb INK notitle, \
         track u 9:($13 != 0 ? sin($10 * pi / 180) : NaN):($13 > 0 ? OUTWARD : INWARD) \
             w p pt 7 ps 0.55 lc rgb variable notitle
}

# --- Graphique L4 : accord de polarité par fenêtre -------------------------------------
set origin 0.53,0.02
set size 0.47,0.66
set lmargin 9
set rmargin 3
unset colorbox
set xdata time
set timefmt "%Y-%m-%dT%H:%M"
set format x "%b\n%Y"
set xrange ["2020-06-15T00:00":"2023-01-15T00:00"]
set yrange [0:100]
set xtics font ",11" ("Jul\n2020" "2020-07-01T00:00", "Jan\n2021" "2021-01-01T00:00", \
    "Jul\n2021" "2021-07-01T00:00", "Jan\n2022" "2022-01-01T00:00", "Jul\n2022" "2022-07-01T00:00", \
    "Jan\n2023" "2023-01-01T00:00")
set ytics 0,20,100 font ",11"
set xlabel ""
set ylabel "Sector polarity agreement (%)" font ",12"
set title "L4 · Solar Orbiter measured vs PFSS predicted" font ",13" tc rgb INK
set grid lc rgb "#e5e4e0"
set key bottom right font ",11" samplen 2
windowNumber(s) = int(s[2:])
isShown = 'windowNumber(strcol(1)) <= nplot'
isOk = 'strcol(2) eq "ok"'
plot summary u 3:(@isShown && @isOk ? $6 : NaN) w lp lw 2.5 pt 7 ps 1.1 lc rgb BLUE t "PFSS prediction", \
     summary u 3:(@isShown && @isOk ? $7 : NaN) w lp lw 1.5 dt 2 pt 6 lc rgb GREY t "majority-polarity baseline", \
     summary u 3:(@isShown && !(@isOk) ? 4 : NaN) w p pt 9 ps 1.4 lc rgb "#fab219" t "insufficient data (no wind speed)", \
     summary u 3:(windowNumber(strcol(1)) == ncur && @isOk ? $6 : NaN) w p pt 6 ps 3 lw 3 lc rgb INK notitle

unset multiplot
