# Operations manual

How the magneto-sun processing chain is organised, run and traced. The chain turns raw
observations (a GONG synoptic magnetogram and Solar Orbiter in situ data) into a
prediction of the interplanetary magnetic sector polarity and a score against the
measurement, window after window over a whole observation campaign.

## Product levels

| Level | Content | Producer | Location |
|---|---|---|---|
| L1 | Raw inputs as delivered: GONG synoptic map (FITS), Solar Orbiter MAG and SWA-PAS (AMDA ASCII), spacecraft ephemeris (JPL Horizons) | `scripts/fetch_*.sh` | `data/` |
| L2 | Standardised, instrument-independent hourly in situ series: B (RTN), radial wind speed, Carrington position | `insitu_l2` | `products/<window>/L2/` |
| L3 | Science products: spherical-harmonic coefficients (with L-curve), PFSS source-surface map and neutral line, measured vs predicted polarity series and window score | `fit_map`, `pfss_map`, `polarity` | `products/<window>/L3/` |
| L4 | Campaign summary: one line per window, and its figure | `scripts/campaign_summary.sh`, `plots/campaign.gp` | `products/` |

Every product `X` has a provenance record `X.meta.json`: level, creation time, code
version (`git describe`, suffixed `-dirty` when the working tree has uncommitted
changes), SHA-256 of the executable, parameters, SHA-256 of every input and of the
product itself. Each processing step also writes its log to `products/<window>/logs/`.
Example:

```json
{
  "product": "products/w12/L3/polarity_score.txt",
  "level": "L3",
  "created_utc": "2026-09-30T20:49:12Z",
  "software": {"version": "a4a1815-dirty", "executable": "build/Release/polarity", "sha256": "1c49e508…"},
  "parameters": "rss=2.5 min_valid_hours=48",
  "inputs": [
    {"path": "products/w12/L3/fit_coeffs.txt", "sha256": "00caa7e2…"},
    {"path": "products/w12/L2/insitu_hourly.txt", "sha256": "5c1fe8db…"}
  ],
  "sha256": "3693c7e1…"
}
```

## Running a campaign

```bash
make build            # compile (Release) and run the unit tests
make fetch            # download the L1 inputs of every window (needs network)
make -j4 -k campaign  # L2, L3 for every window, then the L4 summary
make status           # per-window status and agreement
make replay           # animation of the source surface over the campaign (GIF, MP4)
```

The full 2020–2022 campaign (33 windows) downloads in about 7 minutes and processes in
about 11 seconds on a laptop (`make -j4`).

- **Incremental.** Each product depends on its inputs, on the parameters
  (`config/pipeline.mk`) and on the compiled code: only stale products are recomputed.
  Changing a parameter reprocesses the L3 products, not the downloads nor the L2 series;
  changing the code reprocesses everything, so that every product is traceable to one
  software version.
- **Restartable.** An interrupted run resumes where it stopped; a failed step leaves no
  partial product (`.DELETE_ON_ERROR`).
- **Fault tolerant.** `make -k` keeps processing the other windows when one fails and
  reports the failure. A window without enough valid data is not a failure: its score is
  produced with `status insufficient_data` (threshold `MIN_VALID_HOURS`).
- **Network separated from compute.** Only `make fetch` needs network access. On a
  cluster (ROMEO), run it on the login node, then run `make campaign` on compute nodes.

## Configuration

- `config/pipeline.mk`: processing parameters (spherical-harmonic degree, fitted latitude
  range, source-surface radius, minimum valid hours), version-controlled.
- `config/campaign_*.txt`: the campaign table, one window per line (id, start, end,
  hours, GONG map). Generated once by `scripts/make_campaign.sh`, which picks the first
  GONG map of each window's middle day, then version-controlled:

```bash
scripts/make_campaign.sh 2020-07-14 33 27 > config/campaign_solo_2020_2022.txt
make CAMPAIGN=config/campaign_solo_2020_2022.txt campaign
```

## Replay

`make replay` (`scripts/make_replay.sh`, `plots/replay_frame.gp`) turns the L3
source-surface products of a processed campaign into an animation, one map per window:
radial field at 2.5 solar radii and the neutral line (base of the heliospheric current
sheet). Output: `figures/replay.gif` and `figures/replay.mp4`.

## Continuous integration

`.github/workflows/ci.yml` builds the code on every push and pull request (Debug with
array-bounds checking, and Release), runs all unit tests (Fortran and C++), and dry-runs
the campaign to validate the Makefile and the campaign table.
