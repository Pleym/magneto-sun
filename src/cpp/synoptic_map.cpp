// Lecture des cartes synoptiques horaires GONG et bilan de flux.
#include "synoptic_map.hpp"

#include <fitsio.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
// Écart relatif toléré sur l'étendue d'un axe : CDELT2 est arrondi dans les en-têtes,
// 0.0111111 pour GONG (180 * CDELT2 = 1.999998) et 0.001389 pour HMI (1440 * CDELT2 = 2.00016).
constexpr double SPAN_TOLERANCE = 2e-4;
// Écart toléré sur la position d'un centre de pixel en longitude, en fraction de pixel.
constexpr double PIXEL_TOLERANCE = 1e-3;
// Idem en latitude : l'arrondi de CDELT2 de HMI décale le bord de la grille de 0,06 pixel.
constexpr double LAT_PIXEL_TOLERANCE = 0.1;
constexpr int NO_MISSING_VALUES = 0;

struct FitsCloser {
    void operator()(fitsfile* fptr) const {
        int status = 0;
        fits_close_file(fptr, &status);
    }
};
using FitsPtr = std::unique_ptr<fitsfile, FitsCloser>;

[[noreturn]] void fail(const std::string& path, const std::string& reason) {
    throw std::runtime_error(path + " : " + reason);
}

void checkStatus(int status, const std::string& path, const std::string& context) {
    if (status == 0) {
        return;
    }
    char text[FLEN_STATUS];
    fits_get_errstatus(status, text);
    fail(path, context + " (cfitsio : " + text + ")");
}

FitsPtr openFits(const std::string& path) {
    int status = 0;
    fitsfile* fptr = nullptr;
    fits_open_file(&fptr, path.c_str(), READONLY, &status);
    checkStatus(status, path, "ouverture impossible");
    return FitsPtr(fptr);
}

std::string readString(fitsfile* f, const std::string& path, const char* key) {
    int status = 0;
    char value[FLEN_VALUE];
    fits_read_key(f, TSTRING, key, value, nullptr, &status);
    checkStatus(status, path, std::string("mot-clé ") + key + " absent ou illisible");
    return value;
}

double readDouble(fitsfile* f, const std::string& path, const char* key) {
    int status = 0;
    double value = 0.0;
    fits_read_key(f, TDOUBLE, key, &value, nullptr, &status);
    checkStatus(status, path, std::string("mot-clé ") + key + " absent ou illisible");
    return value;
}

int readInt(fitsfile* f, const std::string& path, const char* key) {
    int status = 0;
    int value = 0;
    fits_read_key(f, TINT, key, &value, nullptr, &status);
    checkStatus(status, path, std::string("mot-clé ") + key + " absent ou illisible");
    return value;
}

bool hasKey(fitsfile* f, const char* key) {
    int status = 0;
    char card[FLEN_CARD];
    fits_read_card(f, key, card, &status);
    fits_clear_errmsg();
    return status == 0;
}

double readDoubleOr(fitsfile* f, const std::string& path, const char* key, double fallback) {
    return hasKey(f, key) ? readDouble(f, path, key) : fallback;
}

int readIntOr(fitsfile* f, const std::string& path, const char* key, int fallback) {
    return hasKey(f, key) ? readInt(f, path, key) : fallback;
}

void requireEqual(const std::string& path, const char* key, const std::string& value,
                  const std::string& expected) {
    if (value != expected) {
        fail(path, std::string(key) + " = '" + value + "', '" + expected + "' attendu");
    }
}

double wrap360(double deg) { return std::fmod(std::fmod(deg, 360.0) + 360.0, 360.0); }

std::pair<int, int> readImageSize(fitsfile* f, const std::string& path) {
    int status = 0;
    int nAxis = 0;
    fits_get_img_dim(f, &nAxis, &status);
    checkStatus(status, path, "dimensions illisibles");
    if (nAxis != 2) {
        fail(path, "image à " + std::to_string(nAxis) + " axe(s), 2 attendus");
    }
    long nAxes[2] = {0, 0};
    fits_get_img_size(f, 2, nAxes, &status);
    checkStatus(status, path, "dimensions illisibles");
    if (nAxes[0] < 1 || nAxes[1] < 1) {
        fail(path, "image vide");
    }
    return {static_cast<int>(nAxes[0]), static_cast<int>(nAxes[1])};
}

struct LongitudeAxis {
    double firstLonDeg;   // premier pixel du fichier, ramené dans [0, 360[
    double stepDeg;       // CDELT1 : négatif si les longitudes décroissent (HMI)
    double offsetPixels;  // centre de la colonne canonique 0, en fraction de pixel
};

// Les cartes commencent à une longitude quelconque (GONG : LONG0 ; HMI : longitudes
// décroissantes, CRVAL1 en « temps de Carrington ») : on relève où tombe la première
// colonne pour les remettre à partir de 0°, dans l'ordre croissant.
LongitudeAxis checkLongitudeAxis(fitsfile* f, const std::string& path, int nLon) {
    requireEqual(path, "CTYPE1", readString(f, path, "CTYPE1"), "CRLN-CEA");
    const double step = readDouble(f, path, "CDELT1");
    if (std::abs(std::abs(step) * nLon / 360.0 - 1.0) > SPAN_TOLERANCE) {
        fail(path, "|CDELT1| * NAXIS1 = " + std::to_string(std::abs(step) * nLon) +
                       "°, 360° attendus");
    }
    const double firstLon = wrap360(readDouble(f, path, "CRVAL1") +
                                    (1.0 - readDouble(f, path, "CRPIX1")) * step);
    const double position = firstLon * nLon / 360.0;
    const double offset = position - std::floor(position + PIXEL_TOLERANCE);
    return {firstLon, step, std::abs(offset) < PIXEL_TOLERANCE ? 0.0 : offset};
}

void checkSineLatitudeAxis(fitsfile* f, const std::string& path, int nLat) {
    requireEqual(path, "CTYPE2", readString(f, path, "CTYPE2"), "CRLT-CEA");
    // PV2_1 est le paramètre lambda de la projection CEA : 1 pour y = sin(latitude).
    const double lambda = readDoubleOr(f, path, "PV2_1", 1.0);
    if (lambda != 1.0) {
        fail(path, "PV2_1 = " + std::to_string(lambda) + ", 1 attendu (sinus de latitude)");
    }
    const double step = readDouble(f, path, "CDELT2");
    if (std::abs(step * nLat / 2.0 - 1.0) > SPAN_TOLERANCE) {
        fail(path, "CDELT2 * NAXIS2 = " + std::to_string(step * nLat) +
                       ", 2 attendu (sinus de latitude de -1 à 1)");
    }
    const double firstSinLat =
        readDouble(f, path, "CRVAL2") + (1.0 - readDouble(f, path, "CRPIX2")) * step;
    const double expected = -1.0 + 1.0 / nLat;
    if (std::abs(firstSinLat - expected) > LAT_PIXEL_TOLERANCE * 2.0 / nLat) {
        fail(path, "première ligne à sin(lat) = " + std::to_string(firstSinLat) + ", " +
                       std::to_string(expected) + " attendu (bord du pôle sud)");
    }
}

std::vector<double> readPixels(fitsfile* f, const std::string& path, int nLon, int nLat) {
    std::vector<double> pixels(static_cast<size_t>(nLon) * nLat);
    long firstPixel[2] = {1, 1};
    int status = 0;
    int hasNull = 0;
    fits_read_pix(f, TDOUBLE, firstPixel, pixels.size(), nullptr, pixels.data(), &hasNull,
                  &status);
    checkStatus(status, path, "lecture des pixels impossible");
    // Pixels manquants acceptés seulement s'ils sont déclarés (MISSVALS, cartes HMI)
    const auto badCount =
        std::count_if(pixels.begin(), pixels.end(), [](double v) { return !std::isfinite(v); });
    const int declared = readIntOr(f, path, "MISSVALS", NO_MISSING_VALUES);
    if (badCount > 0 && badCount != declared) {
        fail(path, std::to_string(badCount) + " pixel(s) non fini(s) (NaN ou infini), " +
                       std::to_string(declared) + " déclaré(s) par MISSVALS");
    }
    return pixels;
}

// Colonnes remises dans l'ordre des longitudes croissantes, la colonne k centrée en
// (k + offset) * 360 / nLon.
std::vector<double> toCanonicalColumns(const std::vector<double>& pixels, int nLon, int nLat,
                                       const LongitudeAxis& axis, const std::string& path) {
    std::vector<int> column(nLon);
    std::vector<bool> isUsed(nLon, false);
    for (int i = 0; i < nLon; ++i) {
        const double position = wrap360(axis.firstLonDeg + i * axis.stepDeg) * nLon / 360.0;
        const int k = static_cast<int>(std::lround(position - axis.offsetPixels) % nLon);
        if (isUsed[k]) {
            fail(path, "colonnes de longitude non alignées sur une grille régulière");
        }
        isUsed[k] = true;
        column[i] = k;
    }
    std::vector<double> canonical(pixels.size());
    for (int j = 0; j < nLat; ++j) {
        const size_t row = static_cast<size_t>(j) * nLon;
        for (int i = 0; i < nLon; ++i) {
            canonical[row + column[i]] = pixels[row + i];
        }
    }
    return canonical;
}

// UT de la carte : MAPDATE/MAPTIME (GONG) ou T_OBS (HMI).
std::string readObservationTime(fitsfile* f, const std::string& path) {
    if (hasKey(f, "MAPDATE")) {
        return readString(f, path, "MAPDATE") + "T" + readString(f, path, "MAPTIME");
    }
    return readString(f, path, "T_OBS");
}

}  // namespace

double SynopticMap::longitudeDeg(int iLon) const {
    return (iLon + lonOffsetPixels) * 360.0 / nLon;
}

double SynopticMap::sinLatitude(int iLat) const { return -1.0 + (iLat + 0.5) * 2.0 / nLat; }

SynopticMap readSynopticMap(const std::string& path) {
    const FitsPtr file = openFits(path);
    fitsfile* f = file.get();
    const auto [nLon, nLat] = readImageSize(f, path);
    const LongitudeAxis axis = checkLongitudeAxis(f, path, nLon);
    checkSineLatitudeAxis(f, path, nLat);
    // 1 Mx/cm^2 = 1 G : les deux unités désignent la même densité de flux
    const std::string unit = readString(f, path, "BUNIT");
    if (unit != "Gauss" && unit != "Mx/cm^2") {
        fail(path, "BUNIT = '" + unit + "', 'Gauss' ou 'Mx/cm^2' attendu");
    }

    return SynopticMap{
        nLon,
        nLat,
        toCanonicalColumns(readPixels(f, path, nLon, nLat), nLon, nLat, axis, path),
        readInt(f, path, "CAR_ROT"),
        readObservationTime(f, path),
        axis.offsetPixels,
    };
}

FluxBalance computeFluxBalance(const SynopticMap& map) {
    // Pixels d'aire égale : chacun couvre 4 pi R^2 / (nLon * nLat). Pixels manquants ignorés.
    const double pixelAreaCm2 =
        4.0 * PI * SOLAR_RADIUS_CM * SOLAR_RADIUS_CM / (static_cast<double>(map.nLon) * map.nLat);
    double net = 0.0;
    double unsignedSum = 0.0;
    for (const double b : map.br) {
        if (!std::isfinite(b)) {
            continue;
        }
        net += b;
        unsignedSum += std::abs(b);
    }
    return {net * pixelAreaCm2, unsignedSum * pixelAreaCm2};
}

void writeRawBinary(const std::vector<double>& values, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(double)));
    if (!out) {
        throw std::runtime_error(path + " : écriture impossible");
    }
}
