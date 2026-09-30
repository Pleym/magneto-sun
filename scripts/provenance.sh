#!/usr/bin/env bash
# Écrit la fiche de traçabilité <produit>.meta.json d'un produit de la chaîne : niveau,
# date de création, version du code (git ; « -dirty » si des modifications ne sont pas
# commitées), exécutable et son empreinte, paramètres, entrées et leurs empreintes,
# empreinte du produit lui-même.
# Usage : scripts/provenance.sh <produit> <niveau> <exécutable> "<paramètres>" <entrées...>
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "usage : $0 <produit> <niveau> <exécutable> \"<paramètres>\" <entrées...>" >&2
  exit 2
fi
product=$1 level=$2 executable=$3 parameters=$4
shift 4

sha256() {
  if command -v sha256sum >/dev/null; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}
version=$(git -C "$(dirname "$0")/.." describe --always --dirty 2>/dev/null || echo unknown)

# ponytail: JSON écrit à la main ; chemins et paramètres sans guillemets ni antislash
{
  printf '{\n'
  printf '  "product": "%s",\n' "$product"
  printf '  "level": "%s",\n' "$level"
  printf '  "created_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "software": {"version": "%s", "executable": "%s", "sha256": "%s"},\n' \
    "$version" "$executable" "$(sha256 "$executable")"
  printf '  "parameters": "%s",\n' "$parameters"
  printf '  "inputs": [\n'
  i=0
  for input in "$@"; do
    i=$((i + 1))
    separator=$([[ $i -lt $# ]] && echo "," || true)
    printf '    {"path": "%s", "sha256": "%s"}%s\n' "$input" "$(sha256 "$input")" "$separator"
  done
  printf '  ],\n'
  printf '  "sha256": "%s"\n' "$(sha256 "$product")"
  printf '}\n'
} >"$product.meta.json"
