#!/usr/bin/env bash
# Télécharge la carte synoptique HMI (SDO) d'une rotation de Carrington dans data/
# (3600 x 1440 pixels, ~21 Mo). À lancer sur la frontale ou en local.
# Usage : scripts/fetch_hmi.sh 2233
set -euo pipefail

cr="${1:?usage : $0 <numéro de rotation de Carrington>}"
if [[ ! "$cr" =~ ^[0-9]{4}$ ]]; then
  echo "numéro de rotation invalide : $cr" >&2
  exit 1
fi
file="hmi.Synoptic_Mr.$cr.fits"
dest="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$dest"
if [[ -f "$dest/$file" ]]; then
  echo "déjà présent : data/$file"
  exit 0
fi
curl -fSL --retry 3 -o "$dest/$file.part" "http://jsoc.stanford.edu/data/hmi/synoptic/$file"
mv "$dest/$file.part" "$dest/$file"
echo "data/$file"
