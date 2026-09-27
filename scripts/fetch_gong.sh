#!/usr/bin/env bash
# Télécharge une carte synoptique horaire GONG (produit zqs) dans data/.
# À lancer sur la frontale ou en local : les nœuds de calcul ne lisent que data/.
# Usage : scripts/fetch_gong.sh mrzqs200701t0014c2232_129.fits.gz
set -euo pipefail

file="${1:?usage : $0 mrzqsAAMMJJtHHMMcRRRR_LLL.fits.gz}"
if [[ ! "$file" =~ ^mrzqs([0-9]{2})([0-9]{2})([0-9]{2})t[0-9]{4}c[0-9]{4}_[0-9]{3}\.fits\.gz$ ]]; then
  echo "nom de carte GONG zqs invalide : $file" >&2
  exit 1
fi
yy=${BASH_REMATCH[1]} mm=${BASH_REMATCH[2]} dd=${BASH_REMATCH[3]}
url="https://gong2.nso.edu/oQR/zqs/20$yy$mm/mrzqs$yy$mm$dd/$file"

dest="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$dest"
if [[ -f "$dest/$file" ]]; then
  echo "déjà présent : data/$file"
  exit 0
fi
curl -fSL --retry 3 -o "$dest/$file.part" "$url"
mv "$dest/$file.part" "$dest/$file"
echo "data/$file"
