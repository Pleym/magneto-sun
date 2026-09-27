// Tests de la Phase 4 : dates, lecture Horizons, interpolation, moyennes horaires,
// spirale de Parker, projection balistique, secteurs.
#include "insitu.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool isOk, const std::string& what) {
    if (!isOk) {
        std::cerr << "ÉCHEC : " << what << '\n';
        ++failures;
    }
}

bool throwsRuntimeError(const std::function<void()>& action) {
    try {
        action();
    } catch (const std::runtime_error&) {
        return true;
    }
    return false;
}

void testIsoTimes() {
    check(parseIsoTime("1970-01-01T00:00:00") == 0.0, "époque Unix");
    check(parseIsoTime("2020-07-01T00:14:00") == 1593562440.0, "2020-07-01T00:14");
    check(std::abs(parseIsoTime("2020-06-20T00:00:29.999") - 1592611229.999) < 1e-6,
          "millisecondes (format AMDA)");
    check(formatIsoTime(1593562440.0) == "2020-07-01T00:14", "formatage");
    check(throwsRuntimeError([] { parseIsoTime("2020-13-01T00:00:00"); }), "mois 13 refusé");
}

void testReadsHorizonsTable() {
    const std::string path =
        (std::filesystem::temp_directory_path() / "magnetosun_horizons.txt").string();
    // Extrait réel : Soleil vu de Solar Orbiter (-144), 2020-07-14
    std::ofstream(path) << "Target body name: Sun (10)\n$$SOE\n"
                           "2459044.500000000, , , 127.492272,   2.225129,  0.62951743566459, 11.0229885,\n"
                           "2459044.541666667, , , 126.986116,   2.215763,  0.62978274839876, 11.0271398,\n"
                           "$$EOE\n";

    const Ephemeris eph = readHorizons(path);
    std::remove(path.c_str());

    check(eph.t.size() == 2, "deux lignes lues");
    check(formatIsoTime(eph.t[0]) == "2020-07-14T00:00", "JD 2459044.5 = 2020-07-14T00:00");
    check(eph.carrLonDeg[0] == 127.492272 && eph.rAu[1] == 0.62978274839876, "colonnes");
}

void testInterpolatesAcrossLongitudeWrap() {
    const Ephemeris eph{{0.0, 3600.0}, {1.0, 359.0}, {2.0, 4.0}, {0.5, 0.7}};

    const Position p = interpolatePosition(eph, 1800.0);

    check(std::abs(p.carrLonDeg) < 1e-9 || std::abs(p.carrLonDeg - 360.0) < 1e-9,
          "1° -> 359° : milieu à 0°, pas à 180°");
    check(std::abs(p.carrLatDeg - 3.0) < 1e-12 && std::abs(p.rAu - 0.6) < 1e-12,
          "latitude et distance interpolées");
    check(throwsRuntimeError([&] { interpolatePosition(eph, 4000.0); }), "hors éphéméride refusé");
}

void testHourlyMeansSkipNaN() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const VectorSeries s{{600, 1200, 1800, 2400},
                         {{{1, 0, 0}}, {{2, 0, 0}}, {{3, 0, 0}}, {{nan, 0, 0}}}};

    const auto means = hourlyMeans(s, 0.0, 2);

    check(means[0][0] == 2.0, "moyenne horaire des valeurs finies");
    check(std::isnan(means[1][0]), "heure sans donnée : NaN");
}

void testParkerSpiralProjection() {
    // À 1 UA et 400 km/s : tan(psi) = Omega r / v
    const double omega = CARRINGTON_RATE_DEG_PER_DAY * 3.14159265358979323846 / 180.0 / 86400.0;
    const double psi = std::atan(omega * AU_KM / 400.0);

    check(std::abs(parkerAlignedField(std::cos(psi), -std::sin(psi), 1.0, 400.0) - 1) < 1e-12,
          "champ le long de la spirale sortante : +|B|");
    check(std::abs(parkerAlignedField(-std::cos(psi), std::sin(psi), 1.0, 400.0) + 1) < 1e-12,
          "spirale entrante : -|B|");
    check(std::abs(parkerAlignedField(std::sin(psi), std::cos(psi), 1.0, 400.0)) < 1e-12,
          "champ perpendiculaire à la spirale : 0");
}

void testBallisticMapping() {
    // 1 UA à 400 km/s : tau = 4,278 jours, soit 60,7° de rotation de Carrington
    const double lon = sourceSurfaceLongitude(10.0, 1.0, 400.0, 2.5);
    check(std::abs(lon - 70.685) < 0.01, "décalage de 60,7° à 1 UA pour 400 km/s");
    check(std::abs(sourceSurfaceLongitude(350.0, 1.0, 400.0, 2.5) - 50.685) < 0.01,
          "longitude ramenée dans [0, 360[");
}

void testSectorIgnoresSwitchbacks() {
    // 48 h sortantes, switchback de 3 h, 45 h sortantes, puis 72 h entrantes
    std::vector<int> hourly(48, 1);
    hourly.insert(hourly.end(), 3, -1);
    hourly.insert(hourly.end(), 45, 1);
    hourly.insert(hourly.end(), 72, -1);

    const std::vector<size_t> changes = polarityChanges(sectorPolarity(hourly, 13));

    check(changes.size() == 1 && changes[0] == 96,
          "un seul changement de secteur, à l'heure 96 (switchback ignoré)");
    check(polarityChanges(hourly).size() == 3, "sans filtrage : 3 changements");
}

}  // namespace

int main() {
    testIsoTimes();
    testReadsHorizonsTable();
    testInterpolatesAcrossLongitudeWrap();
    testHourlyMeansSkipNaN();
    testParkerSpiralProjection();
    testBallisticMapping();
    testSectorIgnoresSwitchbacks();
    if (failures > 0) {
        std::cerr << failures << " vérification(s) en échec\n";
        return 1;
    }
    std::puts("insitu OK");
    return 0;
}
