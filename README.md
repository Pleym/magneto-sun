# magneto-sun

[![CI](https://github.com/Pleym/magneto-sun/actions/workflows/ci.yml/badge.svg)](https://github.com/Pleym/magneto-sun/actions/workflows/ci.yml)

Can a single map of the Sun's surface magnetic field predict the magnetic polarity of the
solar wind measured by a spacecraft? magneto-sun answers with a data processing chain
written in Fortran and C++, from raw observations to scored science products.

![Solar wind source surface, 2020–2022](figures/replay.gif)

*The Sun's magnetic field where the solar wind starts (2.5 solar radii), one map per
27-day window from July 2020 to December 2022, computed by the chain from GONG
magnetograms. Red: field pointing away from the Sun; blue: toward the Sun. Black: the
current sheet between them, the boundary Solar Orbiter crosses when the measured polarity
flips. Nearly flat at solar minimum in 2020, it warps strongly as the Sun becomes more
active. Video version: [figures/replay.mp4](figures/replay.mp4).*

## Results

- **Polarity prediction.** Over 27 windows, the model predicts the measured polarity
  **85 % of the time** on average, against 63 % for a trivial "always the majority
  polarity" baseline. It beats the baseline in 24 windows out of 27.
- **Missing data handled.** 6 windows (Nov 2020 – Apr 2021) have no solar wind speed
  measurement: the chain flags them as `insufficient_data` instead of failing.
- **Validated physics.** The coronal model agrees with the independent solver
  [sunkit-magex](https://github.com/sunpy/sunkit-magex) on 99.9 % of the source surface
  (correlation 1.0000).
- **Performance.** An exact block decomposition of the spherical-harmonic fit gives the
  same solution about **6,000× faster** than the direct least-squares solver (24 s → 4 ms
  at degree 60). A full SDO/HMI map (5.1 million pixels, 520,000 unknowns) is fitted in
  10 s on a laptop.

![Agreement per window, and Solar Orbiter latitude and distance](figures/campaign_summary.png)

## How it is built: a miniature ground segment

Space missions process their data in a *ground segment*: data arrive, go through
successive processing levels, and become science products that scientists can trust.
This project reproduces that organisation at small scale, on real mission data.

- **Processing levels.** L1: raw inputs as delivered (GONG magnetogram, Solar Orbiter
  magnetic field and wind speed from AMDA, spacecraft orbit from JPL Horizons). L2: a
  standardised hourly series, independent of the input formats. L3: science products
  (spherical-harmonic fit, coronal field model, polarity score). L4: campaign summary.
- **Traceability.** Every product comes with a `.meta.json` record: code version,
  checksum of the program, parameters, checksums of all inputs. Any number can be traced
  back to the exact data and code that produced it.
- **Automated operations.** One `make` command runs the whole campaign. Only what changed
  is recomputed, steps run in parallel, an interrupted run resumes where it stopped, and
  a failing window does not stop the others. Downloads are separate from computing, so
  the chain can run on supercomputer nodes without network access.
- **Validation first.** Every algorithm is tested on synthetic data with a known answer
  before real data, then compared with independent references (sunkit-magex, JPL
  Horizons). Build and tests run automatically on every push.
- **Performance engineering.** Profiling-driven optimisation, OpenMP parallelism, and a
  benchmark kit for the ROMEO supercomputer.

Operations manual: [docs/operations.md](docs/operations.md).

## Reproduce

Requirements: CMake, a Fortran and a C++ compiler, cfitsio, LAPACK, gnuplot (and ffmpeg
for the replay).

```bash
make build            # compile and run the tests
make fetch            # download the inputs of the 33 windows (~100 MB)
make -j4 -k campaign  # process every window and build the campaign summary
make replay           # regenerate the replay animation
```
