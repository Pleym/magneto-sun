#include "sh_model.hpp"

#include "magnetosun_fortran.hpp"

#include <cmath>
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

}  // namespace

ShFit fitSynopticMap(const SynopticMap& map, int lmax, double maxAbsLatDeg, double lambda) {
    if (lmax < 0) {
        throw std::runtime_error("lmax doit être positif ou nul");
    }
    const double maxAbsSinLat = std::sin(maxAbsLatDeg * PI / 180.0) + 1e-12;
    std::vector<double> cosTheta;
    std::vector<double> phi;
    std::vector<double> br;
    for (int j = 0; j < map.nLat; ++j) {
        if (std::abs(map.sinLatitude(j)) > maxAbsSinLat) {
            continue;
        }
        for (int i = 0; i < map.nLon; ++i) {
            cosTheta.push_back(map.sinLatitude(j));
            phi.push_back(map.longitudeDeg(i) * PI / 180.0);
            br.push_back(map.br[static_cast<size_t>(j) * map.nLon + i]);
        }
    }

    ShFit fit{zeroCoefficients(lmax), std::vector<double>(N_LAMBDA),
              std::vector<double>(N_LAMBDA), std::vector<double>(N_LAMBDA), 0, 0.0, 0.0,
              static_cast<int>(br.size())};
    int info = 0;
    ms_fit(lmax, fit.nPixels, cosTheta.data(), phi.data(), br.data(), N_LAMBDA, lambda,
           fit.coeffs.g.data(), fit.coeffs.h.data(), fit.lambdas.data(),
           fit.residualNorms.data(), fit.solutionNorms.data(), &fit.cornerIndex, &fit.lambda,
           &fit.conditionNumber, &info);
    if (info == -1) {
        throw std::runtime_error(std::to_string(fit.nPixels) + " pixels pour " +
                                 std::to_string((lmax + 1) * (lmax + 1)) +
                                 " inconnues : système sous-déterminé");
    }
    if (info != 0) {
        throw std::runtime_error("échec LAPACK dans l'ajustement (info = " +
                                 std::to_string(info) + ")");
    }
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
