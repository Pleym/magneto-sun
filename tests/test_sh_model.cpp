// Tests de la frontière C++ / Fortran : rangement des coefficients, ajustement d'une
// carte synthétique, lecture/écriture des coefficients, erreurs.
#include "sh_model.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
int failures = 0;

void check(bool isOk, const std::string& what) {
    if (!isOk) {
        std::cerr << "ÉCHEC : " << what << '\n';
        ++failures;
    }
}

std::string throwMessage(const std::function<void()>& action) {
    try {
        action();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "";
}

double coef(const std::vector<double>& c, int lmax, int l, int m) {
    return c[static_cast<size_t>(m) * (lmax + 1) + l];
}

// B_r = 0.1 + 2 P(1,0) + 0.3 P(1,1) sin(phi) + 0.5 P(2,1) cos(phi), en Schmidt.
SynopticMap knownMap(int nLon, int nLat) {
    SynopticMap map{nLon, nLat, std::vector<double>(static_cast<size_t>(nLon) * nLat), 0, ""};
    for (int j = 0; j < nLat; ++j) {
        const double x = map.sinLatitude(j);
        const double s = std::sqrt(1 - x * x);
        for (int i = 0; i < nLon; ++i) {
            const double phi = map.longitudeDeg(i) * PI / 180.0;
            map.br[static_cast<size_t>(j) * nLon + i] =
                0.1 + 2 * x + 0.3 * s * std::sin(phi) + 0.5 * std::sqrt(3.0) * x * s * std::cos(phi);
        }
    }
    return map;
}

void testFitRecoversKnownCoefficientsInFortranLayout() {
    const SynopticMap map = knownMap(72, 36);

    const ShFit fit = fitSynopticMap(map, 3, 90.0, 0.0);

    const auto& g = fit.coeffs.g;
    const auto& h = fit.coeffs.h;
    check(std::abs(coef(g, 3, 0, 0) - 0.1) < 1e-10, "g(0,0)");
    check(std::abs(coef(g, 3, 1, 0) - 2.0) < 1e-10, "g(1,0)");
    check(std::abs(coef(h, 3, 1, 1) - 0.3) < 1e-10, "h(1,1)");
    check(std::abs(coef(g, 3, 2, 1) - 0.5) < 1e-10, "g(2,1)");
    check(std::abs(coef(g, 3, 3, 2)) < 1e-10, "g(3,2) nul");
    check(fit.nPixels == 72 * 36, "tous les pixels utilisés sans masque");
}

void testSurfaceModelReproducesMapWithoutMonopole() {
    const SynopticMap map = knownMap(72, 36);
    const ShFit fit = fitSynopticMap(map, 3, 90.0, 0.0);

    const std::vector<double> model = pfssBrOnGrid(fit.coeffs, 2.5, 1.0, 72, 36);

    double maxDiff = 0.0;
    for (size_t k = 0; k < model.size(); ++k) {
        maxDiff = std::max(maxDiff, std::abs(model[k] - (map.br[k] - 0.1)));
    }
    check(maxDiff < 1e-10, "modèle à r = 1 : carte moins le monopôle");
}

void testMaskReducesPixelCount() {
    const ShFit fit = fitSynopticMap(knownMap(72, 36), 3, 60.0, 0.0);
    // Centres en sin(lat) = -1 + (j + 1/2)/18 : |sin(lat)| <= sin(60°) pour j = 2..33
    check(fit.nPixels == 72 * 32, "masque |lat| <= 60° : " + std::to_string(fit.nPixels));
}

void testCoefficientsRoundTrip() {
    const ShFit fit = fitSynopticMap(knownMap(72, 36), 3, 90.0, 0.0);
    const std::string path =
        (std::filesystem::temp_directory_path() / "magnetosun_test_coeffs.txt").string();

    writeShCoefficients(fit.coeffs, path, "test");
    const ShCoefficients back = readShCoefficients(path);
    std::remove(path.c_str());

    check(back.lmax == 3 && back.g == fit.coeffs.g && back.h == fit.coeffs.h,
          "coefficients relus à l'identique");
}

void testRejectsInvalidInputs() {
    const std::string path =
        (std::filesystem::temp_directory_path() / "magnetosun_bad_coeffs.txt").string();
    std::ofstream(path) << "# lmax 2\n3 0 1.0 0.0\n";
    const std::string badLine = throwMessage([&] { readShCoefficients(path); });
    std::remove(path.c_str());
    check(badLine.find("ligne invalide") != std::string::npos,
          "l > lmax refusé (reçu : \"" + badLine + "\")");

    const std::string tooFew = throwMessage([] { fitSynopticMap(knownMap(4, 2), 3, 90.0); });
    check(tooFew.find("sous-déterminé") != std::string::npos,
          "moins de pixels que d'inconnues refusé (reçu : \"" + tooFew + "\")");
}

}  // namespace

int main() {
    testFitRecoversKnownCoefficientsInFortranLayout();
    testSurfaceModelReproducesMapWithoutMonopole();
    testMaskReducesPixelCount();
    testCoefficientsRoundTrip();
    testRejectsInvalidInputs();
    if (failures > 0) {
        std::cerr << failures << " vérification(s) en échec\n";
        return 1;
    }
    std::puts("sh_model OK");
    return 0;
}
