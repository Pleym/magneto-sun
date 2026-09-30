#!/usr/bin/env bash
# Vitrine : anime l'évolution de la surface source (B_r du modèle PFSS et nappe de courant)
# sur toutes les fenêtres d'une campagne déjà traitée (produits de « make campaign »),
# une carte par fenêtre. Sortie : GIF + MP4.
# Usage : scripts/make_replay.sh [table de campagne] [dossier des produits] [sortie sans extension]
set -euo pipefail

campaign=${1:-config/campaign_solo_2020_2022.txt}
products=${2:-products}
out=${3:-figures/replay}
MAPS_PER_SECOND=2
HOLD_SECONDS=2.5
frames=$(mktemp -d)
trap 'rm -rf "$frames"' EXIT

frame=0
for w in $(awk '!/^#/ && NF {print $1}' "$campaign"); do
  map=$(awk -v w="$w" '$1 == w {print $5}' "$campaign")
  ss=$products/$w/L3/pfss_ss.bin
  if [[ ! -f "$ss" ]]; then
    echo "produit absent ($ss) : lancer d'abord « make campaign »" >&2
    exit 1
  fi
  # mrzqsAAMMJJtHHMMcRRRR_LLL.fits.gz : date et rotation de Carrington de la carte GONG
  date="20${map:5:2}-${map:7:2}-${map:9:2}"
  cr=$(grep -oE 'c[0-9]{4}' <<<"$map" | tr -d c)
  frame=$((frame + 1))
  gnuplot -e "out='$frames/$(printf %04d $frame).png'; ss='$ss'; \
    neutral='$products/$w/L3/pfss_neutral.txt'; date='$date'; cr='$cr'" plots/replay_frame.gp
done

# Coupes nettes (un fondu superposait deux nappes et deux dates) ; dernière carte tenue
hold="tpad=stop_mode=clone:stop_duration=$HOLD_SECONDS"
ffmpeg -loglevel error -y -framerate "$MAPS_PER_SECOND" -i "$frames/%04d.png" \
  -vf "$hold,scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128:stats_mode=diff[p];[b][p]paletteuse=dither=bayer:bayer_scale=4:diff_mode=rectangle" \
  "$out.gif"
ffmpeg -loglevel error -y -framerate "$MAPS_PER_SECOND" -i "$frames/%04d.png" \
  -vf "$hold,fps=24" -c:v libx264 -pix_fmt yuv420p -crf 20 "$out.mp4"
ls -lh "$out.gif" "$out.mp4"
