// Lit une carte GONG, affiche ses métadonnées et son bilan de flux, et écrit B_r
// en binaire brut pour gnuplot (float64, longitude la plus rapide, sud -> nord).
// Usage : gong_export <carte GONG .fits[.gz]> <sortie .bin>
#include "synoptic_map.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <iterator>
#include <vector>

namespace {

void printSummary(const SynopticMap& map) {
    std::vector<double> finite;
    std::copy_if(map.br.begin(), map.br.end(), std::back_inserter(finite),
                 [](double b) { return std::isfinite(b); });
    const auto [minIt, maxIt] = std::minmax_element(finite.begin(), finite.end());
    const FluxBalance flux = computeFluxBalance(map);
    std::printf("rotation de Carrington : %d\n", map.carringtonRotation);
    std::printf("date de la carte (UT)  : %s\n", map.observationTime.c_str());
    std::printf("grille                 : %d x %d (longitude x sinus de latitude)\n", map.nLon,
                map.nLat);
    std::printf("pixels manquants       : %zu\n", map.br.size() - finite.size());
    std::printf("B_r min / max          : %.1f / %.1f G\n", *minIt, *maxIt);
    std::printf("flux non signé         : %.3e Mx\n", flux.unsignedFluxMx);
    std::printf("flux net               : %.3e Mx (%.2f %% du flux non signé)\n",
                flux.netFluxMx, 100.0 * flux.netFluxMx / flux.unsignedFluxMx);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage : %s <carte GONG .fits[.gz]> <sortie .bin>\n", argv[0]);
        return 2;
    }
    try {
        const SynopticMap map = readSynopticMap(argv[1]);
        printSummary(map);
        writeRawBinary(map.br, argv[2]);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
