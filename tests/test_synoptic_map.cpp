// Tests du lecteur de cartes GONG, sur des fichiers FITS synthétiques écrits
// avec le même en-tête que les vraies cartes horaires mrzqs.
#include "synoptic_map.hpp"

#include <fitsio.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int N_LON = 360;
constexpr int N_LAT = 180;
constexpr double GONG_FIRST_LON_DEG = 129.5;  // comme mrzqs200701t0014c2232_129
constexpr double GONG_CDELT2 = 0.0111111;     // 2/180 tronqué, comme dans les vrais fichiers
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

struct FakeHeader {
    std::string ctype2 = "CRLT-CEA";
    double firstLonDeg = GONG_FIRST_LON_DEG;
    double lonStepDeg = 360.0 / N_LON;
    std::string bunit = "Gauss";
    int missVals = -1;  // < 0 : mot-clé MISSVALS absent
    std::string tObs;   // non vide : T_OBS (HMI) au lieu de MAPDATE/MAPTIME (GONG)
};

void requireOk(int status) {
    if (status != 0) {
        throw std::logic_error("cfitsio a échoué en écrivant la carte de test");
    }
}

// Écrit une carte au format GONG (par défaut) ou HMI ; brFileOrder suit l'ordre du
// fichier, dont la première colonne est à header.firstLonDeg.
void writeFakeGongMap(const std::string& path, const std::vector<double>& brFileOrder,
                      const FakeHeader& header = {}) {
    int status = 0;
    fitsfile* fptr = nullptr;
    fits_create_file(&fptr, ("!" + path).c_str(), &status);
    long naxes[2] = {N_LON, N_LAT};
    fits_create_img(fptr, FLOAT_IMG, 2, naxes, &status);

    auto putString = [&](const char* key, const std::string& value) {
        fits_update_key(fptr, TSTRING, key, const_cast<char*>(value.c_str()), nullptr, &status);
    };
    auto putDouble = [&](const char* key, double value) {
        fits_update_key(fptr, TDOUBLE, key, &value, nullptr, &status);
    };
    int carRot = 2232;
    fits_update_key(fptr, TINT, "CAR_ROT", &carRot, nullptr, &status);
    if (header.tObs.empty()) {
        putString("MAPDATE", "2020-07-01");
        putString("MAPTIME", "00:14");
    } else {
        putString("T_OBS", header.tObs);
    }
    if (header.missVals >= 0) {
        int missVals = header.missVals;
        fits_update_key(fptr, TINT, "MISSVALS", &missVals, nullptr, &status);
    }
    putString("BUNIT", header.bunit);
    putString("CTYPE1", "CRLN-CEA");
    putString("CTYPE2", header.ctype2);
    putDouble("CRPIX1", N_LON / 2 + 0.5);
    putDouble("CRVAL1", header.firstLonDeg + (N_LON / 2 - 0.5) * header.lonStepDeg);
    putDouble("CDELT1", header.lonStepDeg);
    putDouble("CRPIX2", N_LAT / 2 + 0.5);
    putDouble("CRVAL2", 0.0);
    putDouble("CDELT2", GONG_CDELT2);
    putDouble("PV2_1", 1.0);

    long firstPixel[2] = {1, 1};
    fits_write_pix(fptr, TDOUBLE, firstPixel, N_LON * N_LAT,
                   const_cast<double*>(brFileOrder.data()), &status);
    fits_close_file(fptr, &status);
    requireOk(status);
}

// Champ de test : f(longitude du pixel dans le fichier, sinus de latitude).
std::vector<double> fileOrderMap(const std::function<double(double, double)>& field,
                                 const FakeHeader& header = {}) {
    std::vector<double> br(N_LON * N_LAT);
    for (int j = 0; j < N_LAT; ++j) {
        const double sinLat = -1.0 + (j + 0.5) * 2.0 / N_LAT;
        for (int i = 0; i < N_LON; ++i) {
            const double lon = std::fmod(header.firstLonDeg + i * header.lonStepDeg + 360.0, 360.0);
            br[j * N_LON + i] = field(lon, sinLat);
        }
    }
    return br;
}

const std::string MAP_PATH =
    (std::filesystem::temp_directory_path() / "magnetosun_test_map.fits").string();

void testReadsMetadataAndRollsLongitudeToZero() {
    // Arrange : B_r vaut la longitude du pixel
    writeFakeGongMap(MAP_PATH, fileOrderMap([](double lon, double) { return lon; }));

    // Act
    const SynopticMap map = readSynopticMap(MAP_PATH);

    // Assert
    check(map.nLon == N_LON && map.nLat == N_LAT, "dimensions 360 x 180");
    check(map.carringtonRotation == 2232, "CAR_ROT lu");
    check(map.observationTime == "2020-07-01T00:14", "date de la carte : " + map.observationTime);
    check(map.longitudeDeg(0) == 0.5, "centre de la première colonne à 0,5°");
    bool isRolled = true;
    for (int i = 0; i < N_LON; ++i) {
        isRolled = isRolled && map.br[7 * N_LON + i] == map.longitudeDeg(i);
    }
    check(isRolled, "colonne i centrée en longitude (i + 0,5)° après remise en ordre");
}

void testKeepsSineLatitudeOrder() {
    writeFakeGongMap(MAP_PATH, fileOrderMap([](double, double sinLat) { return sinLat; }));

    const SynopticMap map = readSynopticMap(MAP_PATH);

    bool isOrdered = true;
    for (int j = 0; j < N_LAT; ++j) {
        isOrdered = isOrdered && std::abs(map.br[j * N_LON + 42] - map.sinLatitude(j)) < 1e-6;
    }
    check(isOrdered, "ligne j centrée en sinus de latitude -1 + (j + 0,5) * 2/180");
}

void testDipoleHasZeroNetFluxAndKnownUnsignedFlux() {
    // B_r = sin(latitude) gauss : flux net nul, flux non signé 2 pi R^2
    writeFakeGongMap(MAP_PATH, fileOrderMap([](double, double sinLat) { return sinLat; }));

    const FluxBalance flux = computeFluxBalance(readSynopticMap(MAP_PATH));

    const double expected = 2.0 * PI * SOLAR_RADIUS_CM * SOLAR_RADIUS_CM;
    check(std::abs(flux.unsignedFluxMx / expected - 1.0) < 1e-6, "flux non signé du dipôle");
    check(std::abs(flux.netFluxMx) < 1e-6 * expected, "flux net du dipôle nul");
}

void testReadsHmiLikeMap() {
    // Arrange : comme hmi.Synoptic_Mr, longitudes décroissantes centrées sur des degrés
    // entiers, unité Mx/cm^2, date T_OBS, 3 pixels manquants déclarés par MISSVALS
    const FakeHeader hmi{"CRLT-CEA", 359.0, -1.0, "Mx/cm^2", 3, "2020.07.28_21:10:37_TAI"};
    auto br = fileOrderMap([](double lon, double) { return lon; }, hmi);
    br[0] = br[1] = br[2] = std::numeric_limits<double>::quiet_NaN();
    writeFakeGongMap(MAP_PATH, br, hmi);

    // Act
    const SynopticMap map = readSynopticMap(MAP_PATH);

    // Assert
    check(map.longitudeDeg(0) == 0.0, "HMI : colonne 0 centrée sur 0°");
    bool isOrdered = true;
    for (int i = 0; i < N_LON; ++i) {
        isOrdered = isOrdered && map.br[7 * N_LON + i] == map.longitudeDeg(i);
    }
    check(isOrdered, "HMI : longitudes décroissantes remises dans l'ordre croissant");
    const auto nMissing = std::count_if(map.br.begin(), map.br.end(),
                                        [](double b) { return std::isnan(b); });
    check(nMissing == 3 && std::isnan(map.br[N_LON - 1]), "HMI : pixels manquants conservés");
    check(map.observationTime == "2020.07.28_21:10:37_TAI", "HMI : date T_OBS");
}

void testRejectsMissingFile() {
    const std::string message = throwMessage([] { readSynopticMap("/nonexistent/map.fits"); });
    check(message.find("/nonexistent/map.fits") != std::string::npos,
          "fichier absent : erreur qui nomme le fichier (reçu : \"" + message + "\")");
}

void testRejectsNonSineLatitudeGrid() {
    writeFakeGongMap(MAP_PATH, fileOrderMap([](double, double) { return 0.0; }),
                     FakeHeader{"HPLT-TAN", GONG_FIRST_LON_DEG, 360.0 / N_LON, "Gauss", -1, ""});

    const std::string message = throwMessage([] { readSynopticMap(MAP_PATH); });

    check(message.find("CTYPE2") != std::string::npos,
          "grille non Carrington/sinus de latitude refusée (reçu : \"" + message + "\")");
}

void testRejectsNonFiniteValues() {
    auto br = fileOrderMap([](double, double) { return 1.0; });
    br[1234] = std::numeric_limits<double>::quiet_NaN();
    writeFakeGongMap(MAP_PATH, br);

    const std::string message = throwMessage([] { readSynopticMap(MAP_PATH); });

    check(message.find("non fini") != std::string::npos,
          "pixel NaN refusé (reçu : \"" + message + "\")");
}

}  // namespace

int main() {
    testReadsMetadataAndRollsLongitudeToZero();
    testKeepsSineLatitudeOrder();
    testDipoleHasZeroNetFluxAndKnownUnsignedFlux();
    testReadsHmiLikeMap();
    testRejectsMissingFile();
    testRejectsNonSineLatitudeGrid();
    testRejectsNonFiniteValues();
    std::remove(MAP_PATH.c_str());

    if (failures > 0) {
        std::cerr << failures << " vérification(s) en échec\n";
        return 1;
    }
    std::puts("synoptic_map OK");
    return 0;
}
