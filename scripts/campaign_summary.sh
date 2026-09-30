#!/usr/bin/env bash
# Produit L4 : synthèse d'une campagne, une ligne par fenêtre, à partir des scores L3
# (fichiers <fenêtre>/L3/polarity_score.txt, lignes « clé valeur »).
# Usage : scripts/campaign_summary.sh products/*/L3/polarity_score.txt
set -euo pipefail

keys="status start valid_hours hourly_agreement sector_agreement baseline measured_changes"
keys="$keys predicted_changes r_min_au r_max_au lat_min_deg lat_max_deg"
echo "# window $keys"
for score in "$@"; do
  window=$(basename "$(dirname "$(dirname "$score")")")
  awk -v window="$window" -v keys="$keys" '
    { value[$1] = $2 }
    END {
      n = split(keys, key, " ")
      line = window
      for (i = 1; i <= n; i++) line = line " " ((key[i] in value) ? value[key[i]] : "nan")
      print line
    }' "$score"
done
