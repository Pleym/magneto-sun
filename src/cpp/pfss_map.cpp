// Phase 3 : B_r du modèle PFSS sur la surface source et sa ligne neutre, base de la
// nappe de courant héliosphérique. Écrit :
//   <préfixe>_ss.bin        B_r(Rss) sur la grille 360 x 180 (binaire float64)
//   <préfixe>_neutral.txt   points de la ligne neutre : longitude (°), sin(latitude)
//   <préfixe>_spectrum.txt  l, <B_r^2> du degré l à la surface et en Rss (G^2)
// Usage : pfss_map <coeffs.txt> <Rss en rayons solaires> <préfixe>
#include "sh_model.hpp"
#include "synoptic_map.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int N_LON = 360;
constexpr int N_LAT = 180;
constexpr double PI = 3.14159265358979323846;

struct NeutralPoint {
    double lonDeg;
    double sinLat;
};

// Zéros de B_r entre pixels voisins (interpolation linéaire), en longitude
// (avec raccord 360° -> 0°) et en latitude.
std::vector<NeutralPoint> neutralLine(const SynopticMap& grid, const std::vector<double>& br) {
    std::vector<NeutralPoint> points;
    const double dLon = 360.0 / grid.nLon;
    for (int j = 0; j < grid.nLat; ++j) {
        for (int i = 0; i < grid.nLon; ++i) {
            const double b = br[static_cast<size_t>(j) * grid.nLon + i];
            const double east = br[static_cast<size_t>(j) * grid.nLon + (i + 1) % grid.nLon];
            if (b * east < 0) {
                points.push_back({grid.longitudeDeg(i) + dLon * b / (b - east), grid.sinLatitude(j)});
            }
            if (j + 1 < grid.nLat) {
                const double north = br[static_cast<size_t>(j + 1) * grid.nLon + i];
                if (b * north < 0) {
                    const double t = b / (b - north);
                    points.push_back({grid.longitudeDeg(i),
                                      grid.sinLatitude(j) + t * (2.0 / grid.nLat)});
                }
            }
        }
    }
    return points;
}

void writeNeutralLine(const std::vector<NeutralPoint>& points, const std::string& path) {
    std::ofstream out(path);
    out << "# ligne neutre B_r = 0 : longitude (°) sin(latitude)\n";
    for (const NeutralPoint& p : points) {
        out << std::fmod(p.lonDeg, 360.0) << ' ' << p.sinLat << '\n';
    }
    if (!out) {
        throw std::runtime_error(path + " : écriture impossible");
    }
}

// <B_r^2> sur la sphère du degré l = somme sur m de (g^2 + h^2) / (2l + 1) (Schmidt).
double meanSquare(const ShCoefficients& c, int l) {
    double sum = 0.0;
    for (int m = 0; m <= l; ++m) {
        const size_t k = static_cast<size_t>(m) * (c.lmax + 1) + l;
        sum += c.g[k] * c.g[k] + c.h[k] * c.h[k];
    }
    return sum / (2 * l + 1);
}

void writeSpectrum(const ShCoefficients& surface, const ShCoefficients& sourceSurface,
                   const std::string& path) {
    std::ofstream out(path);
    out << "# l <B_r^2>_l surface (G^2) <B_r^2>_l Rss (G^2)\n";
    for (int l = 1; l <= surface.lmax; ++l) {
        out << l << ' ' << meanSquare(surface, l) << ' ' << meanSquare(sourceSurface, l) << '\n';
    }
    if (!out) {
        throw std::runtime_error(path + " : écriture impossible");
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage : %s <coeffs.txt> <Rss> <préfixe>\n", argv[0]);
        return 2;
    }
    try {
        const std::string prefix = argv[3];
        const double rss = std::stod(argv[2]);
        if (rss <= 1.0) {
            throw std::runtime_error("Rss doit dépasser 1 rayon solaire");
        }
        const ShCoefficients coeffs = readShCoefficients(argv[1]);
        const std::vector<double> brSs = pfssBrOnGrid(coeffs, rss, rss, N_LON, N_LAT);
        const SynopticMap grid{N_LON, N_LAT, {}, 0, ""};
        const std::vector<NeutralPoint> neutral = neutralLine(grid, brSs);

        writeRawBinary(brSs, prefix + "_ss.bin");
        writeNeutralLine(neutral, prefix + "_neutral.txt");
        writeSpectrum(pfssCoefficients(coeffs, rss, 1.0), pfssCoefficients(coeffs, rss, rss),
                      prefix + "_spectrum.txt");

        const auto isPositive = [](double b) { return b > 0; };
        const auto maxAbs = std::max_element(brSs.begin(), brSs.end(), [](double a, double b) {
            return std::abs(a) < std::abs(b);
        });
        double maxNeutralLatDeg = 0.0;
        for (const NeutralPoint& p : neutral) {
            maxNeutralLatDeg = std::max(maxNeutralLatDeg, std::abs(std::asin(p.sinLat)) * 180 / PI);
        }
        std::printf("surface source        : Rss = %.2f R_sun\n", rss);
        std::printf("|B_r| max en Rss      : %.4f G\n", std::abs(*maxAbs));
        std::printf("aire où B_r > 0       : %.1f %%\n",
                    100.0 * std::count_if(brSs.begin(), brSs.end(), isPositive) / brSs.size());
        std::printf("ligne neutre          : %zu points, latitude max %.1f°\n", neutral.size(),
                    maxNeutralLatDeg);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
