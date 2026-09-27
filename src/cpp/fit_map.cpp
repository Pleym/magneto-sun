// Phase 2 : ajuste une carte GONG en harmoniques sphériques (moindres carrés et
// Tikhonov ; par défaut, lambda au coin de la courbe en L) et écrit :
//   <préfixe>_coeffs.txt  coefficients de Schmidt de B_r (G)
//   <préfixe>_lcurve.txt  courbe en L : lambda, |A c - b|, |c|
//   <préfixe>_fit.bin     B_r du modèle à la surface (binaire float64, comme gong_export)
// Usage : fit_map <carte.fits[.gz]> <lmax> <|latitude| max ajustée (°)> <préfixe> [lambda]
#include "sh_model.hpp"
#include "synoptic_map.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double POLAR_CAP_DEG = 60.0;  // calottes où l'on compare les champs polaires

void writeLCurve(const ShFit& fit, const std::string& path) {
    std::ofstream out(path);
    out << "# coin " << fit.cornerIndex << " lambda " << fit.lambdas[fit.cornerIndex]
        << "\n# lambda |A c - b| |c|\n";
    for (size_t k = 0; k < fit.lambdas.size(); ++k) {
        out << fit.lambdas[k] << ' ' << fit.residualNorms[k] << ' ' << fit.solutionNorms[k] << '\n';
    }
    if (!out) {
        throw std::runtime_error(path + " : écriture impossible");
    }
}

// Moyenne et RMS de (a - b) sur les lignes de latitude choisies (pixels d'aire égale).
struct RowStats {
    double meanA;
    double meanB;
    double rmsDiff;
};

RowStats compareRows(const SynopticMap& map, const std::vector<double>& model, double latMinDeg,
                     double latMaxDeg) {
    double sumA = 0.0;
    double sumB = 0.0;
    double sumSq = 0.0;
    int count = 0;
    for (int j = 0; j < map.nLat; ++j) {
        const double latDeg = std::asin(map.sinLatitude(j)) * 180.0 / PI;
        if (latDeg < latMinDeg || latDeg > latMaxDeg) {
            continue;
        }
        for (int i = 0; i < map.nLon; ++i) {
            const size_t k = static_cast<size_t>(j) * map.nLon + i;
            sumA += model[k];
            sumB += map.br[k];
            sumSq += (model[k] - map.br[k]) * (model[k] - map.br[k]);
            ++count;
        }
    }
    return {sumA / count, sumB / count, std::sqrt(sumSq / count)};
}

void printSummary(const SynopticMap& map, const ShFit& fit, const std::vector<double>& model,
                  double maxAbsLatDeg, bool isLCurve) {
    const RowStats fitted = compareRows(map, model, -maxAbsLatDeg, maxAbsLatDeg);
    const RowStats north = compareRows(map, model, POLAR_CAP_DEG, 90.0);
    const RowStats south = compareRows(map, model, -90.0, -POLAR_CAP_DEG);
    const int lmax = fit.coeffs.lmax;
    std::printf("carte                 : CR %d, %s UT\n", map.carringtonRotation,
                map.observationTime.c_str());
    std::printf("inconnues             : %d (lmax = %d)\n", (lmax + 1) * (lmax + 1), lmax);
    std::printf("pixels ajustés        : %d / %d (|lat| <= %.0f°)\n", fit.nPixels,
                map.nLon * map.nLat, maxAbsLatDeg);
    std::printf("conditionnement       : %.3e\n", fit.conditionNumber);
    std::printf("lambda                : %.3e (%s)\n", fit.lambda,
                isLCurve ? "coin de la courbe en L" : "imposé");
    std::printf("RMS résidu ajusté     : %.3f G\n", fitted.rmsDiff);
    std::printf("B_r moyen > %.0f°N     : modèle %+.3f G, GONG %+.3f G\n", POLAR_CAP_DEG,
                north.meanA, north.meanB);
    std::printf("B_r moyen > %.0f°S     : modèle %+.3f G, GONG %+.3f G\n", POLAR_CAP_DEG,
                south.meanA, south.meanB);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5 && argc != 6) {
        std::fprintf(stderr,
                     "usage : %s <carte.fits[.gz]> <lmax> <|latitude| max (°)> <préfixe> "
                     "[lambda]\n",
                     argv[0]);
        return 2;
    }
    try {
        const std::string prefix = argv[4];
        const int lmax = std::stoi(argv[2]);
        const double maxAbsLatDeg = std::stod(argv[3]);
        const double lambda = argc == 6 ? std::stod(argv[5]) : -1.0;

        const SynopticMap map = readSynopticMap(argv[1]);
        const ShFit fit = fitSynopticMap(map, lmax, maxAbsLatDeg, lambda);
        const std::vector<double> model =
            pfssBrOnGrid(fit.coeffs, DEFAULT_SOURCE_SURFACE_RADIUS, 1.0, map.nLon, map.nLat);

        writeShCoefficients(fit.coeffs, prefix + "_coeffs.txt",
                            "CR " + std::to_string(map.carringtonRotation) + " " +
                                map.observationTime + ", |lat| <= " + argv[3] +
                                ", lambda = " + std::to_string(fit.lambda));
        writeLCurve(fit, prefix + "_lcurve.txt");
        writeRawBinary(model, prefix + "_fit.bin");
        printSummary(map, fit, model, maxAbsLatDeg, lambda < 0);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
