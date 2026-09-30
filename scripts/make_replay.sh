#!/usr/bin/env bash
# Vitrine : rejoue une campagne déjà traitée (produits de « make campaign ») sous forme de
# tableau de bord d'opérations. Pour chaque fenêtre, la chaîne L1 -> L4 s'allume étape
# par étape, la carte de la surface source apparaît et le graphique de campagne se
# complète. Toutes les valeurs affichées viennent des produits réels. Sortie : GIF + MP4.
# Usage : scripts/make_replay.sh [table de campagne] [dossier des produits] [sortie sans extension]
set -euo pipefail

campaign=${1:-config/campaign_solo_2020_2022.txt}
products=${2:-products}
out=${3:-figures/replay}
FPS=6
HOLD_SECONDS=3
summary=$products/campaign_summary.txt
if [[ ! -f "$summary" ]]; then
  echo "synthèse absente ($summary) : lancer d'abord « make campaign »" >&2
  exit 1
fi
frames=$(mktemp -d)
trap 'rm -rf "$frames"' EXIT

windows=($(awk '!/^#/ && NF {print $1}' "$campaign"))
total=${#windows[@]}
frame=0
header=""
panel_ss="" panel_neutral="" panel_track="" panel_title=""

# render <état L1..L4> <détail L1..L4> <fenêtres affichées> <fenêtre mise en valeur>
render() {
  frame=$((frame + 1))
  gnuplot -e "out='$frames/$(printf %04d $frame).png'; header='$header'; \
    s1='$1'; s2='$2'; s3='$3'; s4='$4'; d1='$5'; d2='$6'; d3='$7'; d4='$8'; \
    ss='$panel_ss'; neutral='$panel_neutral'; track='$panel_track'; maptitle='$panel_title'; \
    summary='$summary'; nplot=$9; ncur=${10}" plots/replay_frame.gp
}

# Accord moyen des fenêtres « ok » jusqu'à la fenêtre n incluse
mean_agreement() {
  awk -v n="$1" '!/^#/ && $2 == "ok" && substr($1, 2) + 0 <= n {s += $6; c++}
                 END {if (c) printf "%d windows scored, mean %.1f %%", c, s / c; else printf "no window scored yet"}' "$summary"
}

for ((i = 1; i <= total; i++)); do
  w=${windows[i - 1]}
  read -r start stop hours map < <(awk -v w="$w" '$1 == w {print $2, $3, $4, $5}' "$campaign")
  tag="${start:0:10}_${stop:0:10}"
  cr=$(grep -oE 'c[0-9]{4}' <<<"$map" | tr -d c)
  l2=$products/$w/L2/insitu_hourly.txt
  score=$products/$w/L3/polarity_score.txt
  header="Solar Orbiter replay · window $i / $total · ${start:0:10} → ${stop:0:10}"

  d1="CR $cr map · $(grep -vc '^#' "data/solo_mag_rtn_$tag.txt") MAG samples"
  d2="$hours h · $(awk '!/^#/ && $3 != "nan"' "$l2" | wc -l | tr -d ' ') with B · $(awk '!/^#/ && $6 != "nan"' "$l2" | wc -l | tr -d ' ') with speed"
  if [[ $(awk '$1 == "status" {print $2}' "$score") == ok ]]; then
    s3=done
    d3=$(awk '$1 == "sector_agreement" {a = $2} $1 == "baseline" {b = $2} END {printf "agreement %.0f %% (baseline %.0f %%)", a, b}' "$score")
  else
    s3=warn
    d3="no solar wind speed: not scored"
  fi
  d4=$(mean_agreement "$i")

  # 1. les entrées L1 sont arrivées, L2 en cours (la carte affichée reste la précédente)
  render done active wait wait "$d1" "" "" "" $((i - 1)) 0
  # 2. produits L3 : la carte de la fenêtre apparaît
  panel_ss=$products/$w/L3/pfss_ss.bin
  panel_neutral=$products/$w/L3/pfss_neutral.txt
  panel_track=$products/$w/L3/polarity_series.txt
  panel_title="L3 · PFSS source surface (CR $cr) and Solar Orbiter footpoints"
  render done done active wait "$d1" "$d2" "" "" $((i - 1)) 0
  # 3. score de la fenêtre versé à la synthèse L4
  render done done "$s3" active "$d1" "$d2" "$d3" "$d4" "$i" "$i"
done

# Image finale tenue quelques secondes
header="Solar Orbiter 2020–2022 · campaign complete · $total windows"
render done done done done "" "" "" "$(mean_agreement "$total")" "$total" 0
for ((k = 1; k < FPS * HOLD_SECONDS; k++)); do
  frame=$((frame + 1))
  cp "$frames/$(printf %04d $((frame - 1))).png" "$frames/$(printf %04d $frame).png"
done

ffmpeg -loglevel error -y -framerate "$FPS" -i "$frames/%04d.png" \
  -vf "scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle" \
  "$out.gif"
ffmpeg -loglevel error -y -framerate "$FPS" -i "$frames/%04d.png" \
  -c:v libx264 -pix_fmt yuv420p -crf 20 "$out.mp4"
ls -lh "$out.gif" "$out.mp4"
