#include "sh_model.hpp"

#include "magnetosun_fortran.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr int N_LAMBDA = 200;  // points de la courbe en L

size_t coefIndex(int lmax, int l, int m) { return static_cast<size_t>(m) * (lmax + 1) + l; }

ShCoefficients zeroCoefficients(int lmax) {
    const size_t n = static_cast<size_t>(lmax + 1) * (lmax + 1);
    return {lmax, std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)};
}

// Points retenus pour l'ajustement : cos(theta), phi (radians), B_r.
struct Samples {
    std::vector<double> cosTheta;
    std::vector<double> phi;
    std::vector<double> br;
};

bool isRowSelected(const SynopticMap& map, int j, double maxAbsSinLat) {
    return std::abs(map.sinLatitude(j)) <= maxAbsSinLat;
}

// Solveur dense : tous les pixels finis des rangées retenues.
Samples selectPixels(const SynopticMap& map, double maxAbsSinLat) {
    Samples s;
    for (int j = 0; j < map.nLat; ++j) {
        for (int i = 0; isRowSelected(map, j, maxAbsSinLat) && i < map.nLon; ++i) {
            const double b = map.br[static_cast<size_t>(j) * map.nLon + i];
            if (std::isfinite(b)) {
                s.cosTheta.push_back(map.sinLatitude(j));
                s.phi.push_back(map.longitudeDeg(i) * PI / 180.0);
                s.br.push_back(b);
            }
        }
    }
    return s;
}

// Solveur par anneaux : les rangées retenues sans pixel manquant (cosTheta par anneau).
Samples selectRings(const SynopticMap& map, double maxAbsSinLat) {
    Samples s;
    for (int j = 0; j < map.nLat; ++j) {
        const auto first = map.br.begin() + static_cast<std::ptrdiff_t>(j) * map.nLon;
        const bool isComplete =
            std::all_of(first, first + map.nLon, [](double b) { return std::isfinite(b); });
        if (isRowSelected(map, j, maxAbsSinLat) && isComplete) {
            s.cosTheta.push_back(map.sinLatitude(j));
            s.br.insert(s.br.end(), first, first + map.nLon);
        }
    }
    return s;
}

void checkFitInfo(int info, int nPixels, int lmax, int nLon) {
    if (info == -1) {
        throw std::runtime_error(std::to_string(nPixels) + " pixels pour lmax = " +
                                 std::to_string(lmax) + " : système sous-déterminé");
    }
    if (info == -2) {
        throw std::runtime_error("lmax = " + std::to_string(lmax) + " trop grand pour " +
                                 std::to_string(nLon) + " pixels en longitude");
    }
    if (info != 0) {
        throw std::runtime_error("échec LAPACK dans l'ajustement (info = " +
                                 std::to_string(info) + ")");
    }
}

}  // namespace

ShFit fitSynopticMap(const SynopticMap& map, int lmax, double maxAbsLatDeg, double lambda,
                     FitSolver solver) {
    if (lmax < 0) {
        throw std::runtime_error("lmax doit être positif ou nul");
    }
    const double maxAbsSinLat = std::sin(maxAbsLatDeg * PI / 180.0) + 1e-12;
    ShFit fit{zeroCoefficients(lmax), std::vector<double>(N_LAMBDA),
              std::vector<double>(N_LAMBDA), std::vector<double>(N_LAMBDA), 0, 0.0, 0.0, 0, {}};
    int info = 0;
    if (solver == FitSolver::Rings) {
        const Samples rings = selectRings(map, maxAbsSinLat);
        const int nRings = static_cast<int>(rings.cosTheta.size());
        fit.nPixels = nRings * map.nLon;
        ms_fit_rings(lmax, map.nLon, nRings, map.longitudeDeg(0) * PI / 180.0,
                     rings.cosTheta.data(), rings.br.data(), N_LAMBDA, lambda,
                     fit.coeffs.g.data(), fit.coeffs.h.data(), fit.lambdas.data(),
                     fit.residualNorms.data(), fit.solutionNorms.data(), &fit.cornerIndex,
                     &fit.lambda, &fit.conditionNumber, &info, fit.stageSeconds.data());
    } else {
        const Samples pixels = selectPixels(map, maxAbsSinLat);
        fit.nPixels = static_cast<int>(pixels.br.size());
        ms_fit(lmax, fit.nPixels, pixels.cosTheta.data(), pixels.phi.data(), pixels.br.data(),
               N_LAMBDA, lambda, fit.coeffs.g.data(), fit.coeffs.h.data(), fit.lambdas.data(),
               fit.residualNorms.data(), fit.solutionNorms.data(), &fit.cornerIndex, &fit.lambda,
               &fit.conditionNumber, &info, fit.stageSeconds.data());
    }
    checkFitInfo(info, fit.nPixels, lmax, map.nLon);
    return fit;
}

ShCoefficients pfssCoefficients(const ShCoefficients& coeffs, double rss, double r) {
    ShCoefficients scaled = zeroCoefficients(coeffs.lmax);
    ms_pfss_coefs(coeffs.lmax, coeffs.g.data(), coeffs.h.data(), rss, r, scaled.g.data(),
                  scaled.h.data());
    return scaled;
}

std::vector<double> pfssBr(const ShCoefficients& coeffs, double rss, double r,
                           const std::vector<double>& cosTheta, const std::vector<double>& phi) {
    if (cosTheta.size() != phi.size()) {
        throw std::logic_error("pfssBr : cosTheta et phi de tailles différentes");
    }
    std::vector<double> values(cosTheta.size());
    ms_pfss_br(coeffs.lmax, coeffs.g.data(), coeffs.h.data(), rss, r,
               static_cast<int>(values.size()), cosTheta.data(), phi.data(), values.data());
    return values;
}

std::vector<double> pfssBrOnGrid(const ShCoefficients& coeffs, double rss, double r, int nLon,
                                 int nLat) {
    const SynopticMap grid{nLon, nLat, {}, 0, ""};
    std::vector<double> cosTheta;
    std::vector<double> phi;
    for (int j = 0; j < nLat; ++j) {
        for (int i = 0; i < nLon; ++i) {
            cosTheta.push_back(grid.sinLatitude(j));
            phi.push_back(grid.longitudeDeg(i) * PI / 180.0);
        }
    }
    return pfssBr(coeffs, rss, r, cosTheta, phi);
}

void writeShCoefficients(const ShCoefficients& coeffs, const std::string& path,
                         const std::string& comment) {
    std::ofstream out(path);
    out << "# " << comment << "\n# lmax " << coeffs.lmax << "\n# l m g h\n";
    out.precision(17);
    for (int l = 0; l <= coeffs.lmax; ++l) {
        for (int m = 0; m <= l; ++m) {
            const size_t k = coefIndex(coeffs.lmax, l, m);
            out << l << ' ' << m << ' ' << coeffs.g[k] << ' ' << coeffs.h[k] << '\n';
        }
    }
    if (!out) {
        throw std::runtime_error(path + " : écriture impossible");
    }
}

ShCoefficients readShCoefficients(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(path + " : lecture impossible");
    }
    int lmax = -1;
    std::string line;
    while (lmax < 0 && std::getline(in, line)) {
        std::istringstream header(line);
        std::string hash;
        std::string key;
        if (header >> hash >> key && hash == "#" && key == "lmax") {
            header >> lmax;
        }
    }
    if (lmax < 0) {
        throw std::runtime_error(path + " : ligne « # lmax N » absente");
    }
    ShCoefficients coeffs = zeroCoefficients(lmax);
    int l = 0;
    int m = 0;
    double g = 0.0;
    double h = 0.0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream row(line);
        if (!(row >> l >> m >> g >> h) || m < 0 || m > l || l > lmax) {
            throw std::runtime_error(path + " : ligne invalide « " + line + " »");
        }
        coeffs.g[coefIndex(lmax, l, m)] = g;
        coeffs.h[coefIndex(lmax, l, m)] = h;
    }
    return coeffs;
}
