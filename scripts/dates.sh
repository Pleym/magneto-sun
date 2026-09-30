# Conversions de dates UTC portables, à inclure avec « source » : date GNU (Linux,
# ROMEO) ou date BSD (macOS). Format ISO : AAAA-MM-JJTHH:MM:SS.
if date -u -d @0 >/dev/null 2>&1; then
  to_epoch() { date -u -d "${1/T/ }" +%s; }
  to_iso() { date -u -d "@$1" +%Y-%m-%dT%H:%M:%S; }
else
  to_epoch() { date -u -j -f %Y-%m-%dT%H:%M:%S "$1" +%s; }
  to_iso() { date -u -r "$1" +%Y-%m-%dT%H:%M:%S; }
fi
