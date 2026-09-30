// Exploitation, niveau L2 -> L3 : compare la polarité mesurée (série in situ L2) à celle
// que prédit le modèle PFSS (coefficients L3) au pied balistique du vent mesuré. Écrit :
//   <préfixe>_series.txt  série horaire :
//     1 temps (s depuis 1970, UTC)   2 date   3-5 B_R, B_T, B_N (nT)   6 v_R (km/s)
//     7 r (UA)   8 longitude de la sonde (°)   9 longitude du pied (°)   10 latitude (°)
//     11 B le long de la spirale (nT)   12 polarité horaire   13 polarité de secteur
//     14 B_r(Rss) prédit (G)   15 polarité prédite          (NaN / 0 : indéterminé)
//   <préfixe>_score.txt   score de la fenêtre, lignes « clé valeur »
// Usage : polarity <coeffs.txt> <Rss> <série L2> <heures valides minimum> <préfixe>
#include "insitu.hpp"
#include "sh_model.hpp"

#include <algorithm>
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

struct HourAnalysis {
    double footLonDeg = NaN;
    double bParker = NaN;
    int measured = 0;
    int sector = 0;
    double brSourceSurface = NaN;
    int predicted = 0;
};

int signOf(double x) { return std::isfinite(x) ? (x > 0) - (x < 0) : 0; }

std::vector<HourAnalysis> measure(const std::vector<HourlyRecord>& series, double rss) {
    std::vector<HourAnalysis> hours(series.size());
    std::vector<int> measured(series.size(), 0);
    for (size_t k = 0; k < series.size(); ++k) {
        const HourlyRecord& h = series[k];
        if (std::isfinite(h.vR) && std::isfinite(h.bRtn[0])) {
            hours[k].bParker = parkerAlignedField(h.bRtn[0], h.bRtn[1], h.position.rAu, h.vR);
            hours[k].footLonDeg =
                sourceSurfaceLongitude(h.position.carrLonDeg, h.position.rAu, h.vR, rss);
            hours[k].measured = measured[k] = signOf(hours[k].bParker);
        }
    }
    const std::vector<int> sector = sectorPolarity(measured, SECTOR_WINDOW_HOURS);
    for (size_t k = 0; k < hours.size(); ++k) {
        hours[k].sector = hours[k].measured != 0 ? sector[k] : 0;
    }
    return hours;
}

void predict(std::vector<HourAnalysis>& hours, const std::vector<HourlyRecord>& series,
             const ShCoefficients& coeffs, double rss) {
    std::vector<double> cosTheta;
    std::vector<double> phi;
    std::vector<size_t> index;
    for (size_t k = 0; k < hours.size(); ++k) {
        if (std::isfinite(hours[k].footLonDeg)) {
            cosTheta.push_back(std::sin(series[k].position.carrLatDeg * PI / 180.0));
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

FILE* openOutput(const std::string& path) {
    FILE* out = std::fopen(path.c_str(), "w");
    if (out == nullptr) {
        throw std::runtime_error(path + " : écriture impossible");
    }
    return out;
}

void writeSeries(const std::vector<HourlyRecord>& series, const std::vector<HourAnalysis>& hours,
                 const std::string& path) {
    FILE* out = openOutput(path);
    std::fprintf(out, "# t date B_R B_T B_N v_R r lon_sc lon_pied lat B_spirale pol secteur "
                      "Br_Rss pol_predite\n");
    for (size_t k = 0; k < series.size(); ++k) {
        const HourlyRecord& s = series[k];
        const HourAnalysis& h = hours[k];
        std::fprintf(out, "%.0f %s %.3f %.3f %.3f %.1f %.5f %.3f %.3f %.3f %.3f %d %d %.6g %d\n",
                     s.t, formatIsoTime(s.t).c_str(), s.bRtn[0], s.bRtn[1], s.bRtn[2], s.vR,
                     s.position.rAu, s.position.carrLonDeg, h.footLonDeg, s.position.carrLatDeg,
                     h.bParker, h.measured, h.sector, h.brSourceSurface, h.predicted);
    }
    std::fclose(out);
}

// Score et contexte de la fenêtre (distance et latitude de la sonde aux heures mesurées).
void writeScore(const PolarityScore& score, const std::vector<HourlyRecord>& series,
                const std::vector<HourAnalysis>& hours, double rss, const std::string& path) {
    double rMin = NaN, rMax = NaN, latMin = NaN, latMax = NaN;
    for (size_t k = 0; k < series.size(); ++k) {
        if (hours[k].measured == 0) {
            continue;
        }
        const Position& p = series[k].position;
        rMin = std::isnan(rMin) ? p.rAu : std::min(rMin, p.rAu);
        rMax = std::isnan(rMax) ? p.rAu : std::max(rMax, p.rAu);
        latMin = std::isnan(latMin) ? p.carrLatDeg : std::min(latMin, p.carrLatDeg);
        latMax = std::isnan(latMax) ? p.carrLatDeg : std::max(latMax, p.carrLatDeg);
    }
    FILE* out = openOutput(path);
    std::fprintf(out,
                 "status %s\nstart %s\nhours %zu\nrss %.2f\nvalid_hours %d\n"
                 "hourly_agreement %.2f\nsector_agreement %.2f\nbaseline %.2f\n"
                 "measured_changes %zu\npredicted_changes %zu\n"
                 "r_min_au %.4f\nr_max_au %.4f\nlat_min_deg %.3f\nlat_max_deg %.3f\n",
                 score.status.c_str(), formatIsoTime(series.front().t - 1800.0).c_str(),
                 series.size(), rss, score.validHours, score.hourlyAgreement,
                 score.sectorAgreement, score.baseline, score.measuredChanges,
                 score.predictedChanges, rMin, rMax, latMin, latMax);
    std::fclose(out);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr,
                     "usage : %s <coeffs.txt> <Rss> <série L2> <heures valides minimum> "
                     "<préfixe>\n",
                     argv[0]);
        return 2;
    }
    try {
        const double rss = std::stod(argv[2]);
        const int minValidHours = std::stoi(argv[4]);
        const std::string prefix = argv[5];
        if (rss <= 1.0) {
            throw std::runtime_error("Rss doit dépasser 1 rayon solaire");
        }
        const std::vector<HourlyRecord> series = readHourlySeries(argv[3]);
        std::vector<HourAnalysis> hours = measure(series, rss);
        predict(hours, series, readShCoefficients(argv[1]), rss);

        std::vector<int> measured, sector, predicted;
        for (const HourAnalysis& h : hours) {
            measured.push_back(h.measured);
            sector.push_back(h.sector);
            predicted.push_back(h.predicted);
        }
        const PolarityScore score = scorePolarity(measured, sector, predicted, minValidHours);
        writeSeries(series, hours, prefix + "_series.txt");
        writeScore(score, series, hours, rss, prefix + "_score.txt");
        std::printf("statut                    : %s\n", score.status.c_str());
        std::printf("heures valides            : %d / %zu\n", score.validHours, series.size());
        std::printf("accord de secteur         : %.1f %% (référence triviale %.1f %%)\n",
                    score.sectorAgreement, score.baseline);
        std::printf("changements de secteur    : %zu mesurés, %zu prédits\n",
                    score.measuredChanges, score.predictedChanges);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
