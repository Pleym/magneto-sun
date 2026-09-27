#!/usr/bin/env bash
# Télécharge dans data/ les données Solar Orbiter de la Phase 4 :
#   - champ magnétique L2 en RTN, moyennes 1 min (MAG ; AMDA : solo_b_rtn)
#   - vitesse du vent en RTN, moyennes horaires (SWA-PAS ; AMDA : pas_momgr1_v_rtn)
#   - position de la sonde en coordonnées de Carrington, pas de 1 h (JPL Horizons)
# À lancer sur la frontale ou en local : les nœuds de calcul ne lisent que data/.
# Usage : scripts/fetch_solo.sh 2020-07-14T00:00:00 2020-08-10T00:00:00
set -euo pipefail

usage="usage : $0 <début AAAA-MM-JJTHH:MM:SS> <fin AAAA-MM-JJTHH:MM:SS>"
start="${1:?$usage}"
stop="${2:?$usage}"
for t in "$start" "$stop"; do
  if [[ ! "$t" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}$ ]]; then
    echo "date invalide : $t" >&2
    exit 1
  fi
done
dest="$(cd "$(dirname "$0")/.." && pwd)/data"
mkdir -p "$dest"
tag="${start:0:10}_${stop:0:10}"

# amda <paramètre> <fichier de sortie> [échantillonnage en s]
amda() {
  local token response url
  token=$(curl -fsS "https://amda.irap.omp.eu/php/rest/auth.php")
  response=$(curl -fsS "https://amda.irap.omp.eu/php/rest/getParameter.php?startTime=$start&stopTime=$stop&parameterID=$1&token=$token&outputFormat=ASCII&timeFormat=ISO8601${3:+&sampling=$3}")
  url=$(echo "$response" | sed 's#\\/#/#g' | grep -oE 'https://[^"]+' | head -1 || true)
  if [[ -z "$url" ]]; then
    echo "AMDA a refusé la requête $1 : $response" >&2
    exit 1
  fi
  curl -fsS -o "$dest/$2" "$url"
  echo "data/$2"
}

amda solo_b_rtn "solo_mag_rtn_$tag.txt"
# ponytail: moyenne horaire calculée par AMDA ; la vitesse du vent varie lentement
amda pas_momgr1_v_rtn "solo_pas_v_$tag.txt" 3600

# Soleil (10) vu de Solar Orbiter (-144) : le point sous la sonde donne sa longitude
# et sa latitude de Carrington (quantité 14), delta sa distance (quantité 20).
curl -fsS -G "https://ssd.jpl.nasa.gov/api/horizons.api" \
  --data-urlencode "format=text" --data-urlencode "COMMAND='10'" \
  --data-urlencode "CENTER='500@-144'" --data-urlencode "MAKE_EPHEM=YES" \
  --data-urlencode "EPHEM_TYPE=OBSERVER" --data-urlencode "START_TIME='${start/T/ }'" \
  --data-urlencode "STOP_TIME='${stop/T/ }'" --data-urlencode "STEP_SIZE='1h'" \
  --data-urlencode "QUANTITIES='14,20'" --data-urlencode "CAL_FORMAT=JD" \
  --data-urlencode "CSV_FORMAT=YES" -o "$dest/solo_horizons_$tag.txt"
if ! grep -q '^\$\$SOE' "$dest/solo_horizons_$tag.txt"; then
  echo "réponse Horizons sans table, voir data/solo_horizons_$tag.txt" >&2
  exit 1
fi
echo "data/solo_horizons_$tag.txt"
