# magneto-sun

From a solar surface magnetogram to the solar wind measured in situ: can a
potential-field source-surface (PFSS) model, driven by a single GONG synoptic map,
predict the magnetic sector polarity seen by Solar Orbiter?

Fortran (numerical kernels) + C++ (I/O and data processing), built with CMake.

## Results

_Work in progress. This section will present the final results: key figures,
polarity agreement score, and performance and scaling on the ROMEO supercomputer._

## Reproduce

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```
