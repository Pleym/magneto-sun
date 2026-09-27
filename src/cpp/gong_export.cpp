// Lit une carte GONG, affiche ses métadonnées et son bilan de flux, et écrit B_r
// en binaire brut pour gnuplot (float64, longitude la plus rapide, sud -> nord).
// Usage : gong_export <carte GONG .fits[.gz]> <sortie .bin>
#include "synoptic_map.hpp"

#include <algorithm>
#include <cstdio>
#include <exception>

namespace {

void printSummary(const SynopticMap& map) {
    const auto [minIt, maxIt] = std::minmax_element(map.br.begin(), map.br.end());
    const FluxBalance flux = computeFluxBalance(map);
    std::printf("rotation de Carrington : %d\n", map.carringtonRotation);
    std::printf("date de la carte (UT)  : %s\n", map.observationTime.c_str());
    std::printf("grille                 : %d x %d (longitude x sinus de latitude)\n", map.nLon,
                map.nLat);
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
        const SynopticMap map = readGongMap(argv[1]);
        printSummary(map);
        writeRawBinary(map.br, argv[2]);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
