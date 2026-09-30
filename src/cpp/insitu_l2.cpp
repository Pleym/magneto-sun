// Exploitation, niveau L1 -> L2 : fusionne les données in situ brutes (MAG et SWA-PAS
// d'AMDA, éphéméride JPL Horizons) en une série horaire standardisée (HourlyRecord).
// Usage : insitu_l2 <mag.txt> <pas.txt> <horizons.txt> <début ISO> <heures> <sortie L2>
#include "insitu.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 7) {
        std::fprintf(stderr,
                     "usage : %s <mag.txt> <pas.txt> <horizons.txt> <début ISO> <heures> "
                     "<sortie L2>\n",
                     argv[0]);
        return 2;
    }
    try {
        const int nHours = std::stoi(argv[5]);
        if (nHours < 1) {
            throw std::runtime_error("le nombre d'heures doit être positif");
        }
        const std::vector<HourlyRecord> series =
            buildHourlySeries(readAmdaVectors(argv[1]), readAmdaVectors(argv[2]),
                              readHorizons(argv[3]), parseIsoTime(argv[4]), nHours);
        writeHourlySeries(series, argv[6]);

        int withB = 0;
        int withV = 0;
        for (const HourlyRecord& h : series) {
            withB += std::isfinite(h.bRtn[0]);
            withV += std::isfinite(h.vR);
        }
        std::printf("heures                : %zu (%s -> %s)\n", series.size(),
                    formatIsoTime(series.front().t).c_str(), formatIsoTime(series.back().t).c_str());
        std::printf("heures avec B (MAG)   : %d\n", withB);
        std::printf("heures avec v (PAS)   : %d\n", withV);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
