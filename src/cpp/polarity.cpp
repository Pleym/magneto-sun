// Phase 4 : compare la polarité du champ mesurée par Solar Orbiter à celle que prédit
// le modèle PFSS au pied balistique du vent mesuré. Écrit une série horaire :
//   1 temps (s depuis 1970, UTC)   2 date   3-5 B_R, B_T, B_N (nT)   6 v_R (km/s)
//   7 r (UA)   8 longitude de la sonde (°)   9 longitude du pied (°)   10 latitude (°)
//   11 B le long de la spirale (nT)   12 polarité horaire   13 polarité de secteur
//   14 B_r(Rss) prédit (G)   15 polarité prédite          (NaN / 0 : indéterminé)
// Usage : polarity <coeffs.txt> <Rss> <mag.txt> <pas.txt> <horizons.txt> <début ISO>
//                  <nombre d'heures> <sortie.txt>
#include "insitu.hpp"
#include "sh_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
// Fenêtre du vote majoritaire (±6 h) : efface les switchbacks, garde les secteurs.
constexpr int SECTOR_WINDOW_HOURS = 13;

struct HourRecord {
    double t;
    std::array<double, 3> b;
    double vR;
    Position position;
    double footLonDeg = NaN;
    double bParker = NaN;
    int measured = 0;
    int sector = 0;
    double brSourceSurface = NaN;
    int predicted = 0;
};

int signOf(double x) { return std::isfinite(x) ? (x > 0) - (x < 0) : 0; }

std::vector<HourRecord> measure(const VectorSeries& mag, const VectorSeries& pas,
                                const Ephemeris& eph, double t0, int nHours, double rss) {
    const auto b = hourlyMeans(mag, t0, nHours);
    const auto v = hourlyMeans(pas, t0, nHours);
    std::vector<HourRecord> hours;
    std::vector<int> measured;
    for (int k = 0; k < nHours; ++k) {
        const double t = t0 + (k + 0.5) * 3600.0;
        HourRecord h{t, b[k], v[k][0], interpolatePosition(eph, t)};
        if (std::isfinite(h.vR) && std::isfinite(h.b[0])) {
            h.bParker = parkerAlignedField(h.b[0], h.b[1], h.position.rAu, h.vR);
            h.footLonDeg = sourceSurfaceLongitude(h.position.carrLonDeg, h.position.rAu, h.vR, rss);
            h.measured = signOf(h.bParker);
        }
        measured.push_back(h.measured);
        hours.push_back(h);
    }
    const std::vector<int> sector = sectorPolarity(measured, SECTOR_WINDOW_HOURS);
    for (size_t k = 0; k < hours.size(); ++k) {
        hours[k].sector = hours[k].measured != 0 ? sector[k] : 0;
    }
    return hours;
}

void predict(std::vector<HourRecord>& hours, const ShCoefficients& coeffs, double rss) {
    std::vector<double> cosTheta;
    std::vector<double> phi;
    std::vector<size_t> index;
    for (size_t k = 0; k < hours.size(); ++k) {
        if (std::isfinite(hours[k].footLonDeg)) {
            cosTheta.push_back(std::sin(hours[k].position.carrLatDeg * PI / 180.0));
            phi.push_back(hours[k].footLonDeg * PI / 180.0);
            index.push_back(k);
        }
    }
    const std::vector<double> br = pfssBr(coeffs, rss, rss, cosTheta, phi);
    for (size_t n = 0; n < index.size(); ++n) {
        hours[index[n]].brSourceSurface = br[n];
        hours[index[n]].predicted = signOf(br[n]);
    }
}

void writeSeries(const std::vector<HourRecord>& hours, const std::string& path) {
    FILE* out = std::fopen(path.c_str(), "w");
    if (out == nullptr) {
        throw std::runtime_error(path + " : écriture impossible");
    }
    std::fprintf(out, "# t date B_R B_T B_N v_R r lon_sc lon_pied lat B_spirale pol secteur "
                      "Br_Rss pol_predite\n");
    for (const HourRecord& h : hours) {
        std::fprintf(out, "%.0f %s %.3f %.3f %.3f %.1f %.5f %.3f %.3f %.3f %.3f %d %d %.6g %d\n",
                     h.t, formatIsoTime(h.t).c_str(), h.b[0], h.b[1], h.b[2], h.vR,
                     h.position.rAu, h.position.carrLonDeg, h.footLonDeg,
                     h.position.carrLatDeg, h.bParker, h.measured, h.sector,
                     h.brSourceSurface, h.predicted);
    }
    std::fclose(out);
}

void printScores(const std::vector<HourRecord>& hours) {
    int valid = 0, hourlyHits = 0, sectorHits = 0, outward = 0;
    std::vector<int> sector;
    std::vector<int> predicted;
    for (const HourRecord& h : hours) {
        sector.push_back(h.sector);
        predicted.push_back(h.predicted);
        if (h.sector == 0 || h.predicted == 0) {
            continue;
        }
        ++valid;
        hourlyHits += h.measured == h.predicted;
        sectorHits += h.sector == h.predicted;
        outward += h.sector > 0;
    }
    const double baseline = 100.0 * std::max(outward, valid - outward) / valid;
    std::printf("heures valides            : %d / %zu\n", valid, hours.size());
    std::printf("accord horaire            : %.1f %%\n", 100.0 * hourlyHits / valid);
    std::printf("accord de secteur         : %.1f %%\n", 100.0 * sectorHits / valid);
    std::printf("référence triviale        : %.1f %% (toujours la polarité majoritaire)\n",
                baseline);
    std::printf("changements de secteur    :\n");
    for (const size_t k : polarityChanges(sector)) {
        std::printf("  mesuré  %s\n", formatIsoTime(hours[k].t).c_str());
    }
    for (const size_t k : polarityChanges(predicted)) {
        std::printf("  prédit  %s\n", formatIsoTime(hours[k].t).c_str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 9) {
        std::fprintf(stderr,
                     "usage : %s <coeffs.txt> <Rss> <mag.txt> <pas.txt> <horizons.txt> "
                     "<début ISO> <heures> <sortie.txt>\n",
                     argv[0]);
        return 2;
    }
    try {
        const double rss = std::stod(argv[2]);
        const int nHours = std::stoi(argv[7]);
        if (rss <= 1.0 || nHours < 1) {
            throw std::runtime_error("Rss doit dépasser 1 et le nombre d'heures être positif");
        }
        std::vector<HourRecord> hours = measure(readAmdaVectors(argv[3]), readAmdaVectors(argv[4]),
                                                readHorizons(argv[5]), parseIsoTime(argv[6]),
                                                nHours, rss);
        predict(hours, readShCoefficients(argv[1]), rss);
        writeSeries(hours, argv[8]);
        printScores(hours);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
