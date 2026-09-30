#!/usr/bin/env bash
# Nom de la première carte synoptique horaire GONG (zqs) d'une journée ; si la journée
# n'en a pas (panne du réseau GONG), essaie les jours suivants.
# Usage : scripts/find_gong_map.sh 2020-07-27
set -euo pipefail

MAX_DAYS=3
day="${1:?usage : $0 AAAA-MM-JJ}"
if [[ ! "$day" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
  echo "date invalide : $day" >&2
  exit 1
fi
source "$(dirname "$0")/dates.sh"

epoch=$(to_epoch "${day}T00:00:00")
for ((i = 0; i < MAX_DAYS; i++)); do
  d=$(to_iso $((epoch + i * 86400)))
  yy=${d:2:2} mm=${d:5:2} dd=${d:8:2}
  name=$(curl -fsS --retry 3 "https://gong2.nso.edu/oQR/zqs/20$yy$mm/mrzqs$yy$mm$dd/" 2>/dev/null |
    grep -oE "mrzqs$yy$mm${dd}t[0-9]{4}c[0-9]{4}_[0-9]{3}\.fits\.gz" | head -1 || true)
  if [[ -n "$name" ]]; then
    echo "$name"
    exit 0
  fi
done
echo "aucune carte GONG entre le $day et les $((MAX_DAYS - 1)) jours suivants" >&2
exit 1
