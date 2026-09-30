#!/usr/bin/env bash
# Génère la table d'une campagne : fenêtres consécutives de N jours, chacune avec la carte
# GONG du milieu de fenêtre. La table produite est versionnée dans config/ : c'est la
# configuration de la campagne, relue par le Makefile.
# Usage : scripts/make_campaign.sh <début AAAA-MM-JJ> <nombre de fenêtres> <jours par fenêtre>
#         > config/campaign_xxx.txt
set -euo pipefail

usage="usage : $0 <début AAAA-MM-JJ> <nombre de fenêtres> <jours par fenêtre>"
start="${1:?$usage}"
count="${2:?$usage}"
days="${3:?$usage}"
dir="$(dirname "$0")"
source "$dir/dates.sh"

echo "# Campagne magnetosun : $count fenêtres de $days jours à partir du $start"
echo "# id début (UTC) fin (UTC) heures carte GONG du milieu de fenêtre (L1)"
t=$(to_epoch "${start}T00:00:00")
for ((i = 1; i <= count; i++)); do
  begin=$(to_iso "$t")
  end=$(to_iso $((t + days * 86400)))
  middle=$(to_iso $((t + days * 86400 / 2)))
  printf "w%02d %s %s %d %s\n" "$i" "$begin" "$end" $((days * 24)) \
    "$("$dir/find_gong_map.sh" "${middle:0:10}")"
  t=$((t + days * 86400))
done
