#include "insitu.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double SECONDS_PER_DAY = 86400.0;
constexpr double UNIX_EPOCH_JD = 2440587.5;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

double carringtonRateRadPerSecond() {
    return CARRINGTON_RATE_DEG_PER_DAY * PI / 180.0 / SECONDS_PER_DAY;
}

// Jours depuis 1970-01-01 du calendrier grégorien (algorithme de H. Hinnant).
long daysFromCivil(int y, int m, int d) {
    y -= m <= 2;
    const long era = (y >= 0 ? y : y - 399) / 400;
    const int yoe = y - static_cast<int>(era * 400);
    const int doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const int doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + doe - 719468;
}

void civilFromDays(long z, int& y, int& m, int& d) {
    z += 719468;
    const long era = (z >= 0 ? z : z - 146096) / 146097;
    const int doe = static_cast<int>(z - era * 146097);
    const int yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const int doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const int mp = (5 * doy + 2) / 153;
    d = doy - (153 * mp + 2) / 5 + 1;
    m = mp + (mp < 10 ? 3 : -9);
    y = static_cast<int>(yoe + era * 400) + (m <= 2);
}

double wrap360(double deg) { return std::fmod(std::fmod(deg, 360.0) + 360.0, 360.0); }

std::vector<std::string> splitCsv(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) {
        fields.push_back(field);
    }
    return fields;
}

}  // namespace

double parseIsoTime(const std::string& text) {
    int y = 0;
    int mo = 0;
    int d = 0;
    int h = 0;
    int mi = 0;
    double s = 0.0;
    const bool isParsed =
        std::sscanf(text.c_str(), "%4d-%2d-%2dT%2d:%2d:%lf", &y, &mo, &d, &h, &mi, &s) == 6;
    if (!isParsed || mo < 1 || mo > 12 || d < 1 || d > 31 || h > 23 || mi > 59 || s < 0 ||
        s >= 61) {
        throw std::runtime_error("date ISO invalide : « " + text + " »");
    }
    return daysFromCivil(y, mo, d) * SECONDS_PER_DAY + h * 3600.0 + mi * 60.0 + s;
}

std::string formatIsoTime(double t) {
    const long totalMinutes = std::lround(t / 60.0);
    const long days = static_cast<long>(std::floor(totalMinutes / 1440.0));
    const long minuteOfDay = totalMinutes - days * 1440;
    int y = 0;
    int m = 0;
    int d = 0;
    civilFromDays(days, y, m, d);
    char text[32];
    std::snprintf(text, sizeof text, "%04d-%02d-%02dT%02ld:%02ld", y, m, d, minuteOfDay / 60,
                  minuteOfDay % 60);
    return text;
}

VectorSeries readAmdaVectors(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(path + " : lecture impossible");
    }
    VectorSeries series;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream row(line);
        std::string time;
        std::array<std::string, 3> text;
        if (!(row >> time >> text[0] >> text[1] >> text[2])) {
            throw std::runtime_error(path + " : ligne invalide « " + line + " »");
        }
        // std::stod accepte « NaN », valeur de remplissage d'AMDA
        series.t.push_back(parseIsoTime(time));
        series.v.push_back({std::stod(text[0]), std::stod(text[1]), std::stod(text[2])});
    }
    if (series.t.empty()) {
        throw std::runtime_error(path + " : aucune donnée");
    }
    return series;
}

Ephemeris readHorizons(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error(path + " : lecture impossible");
    }
    Ephemeris eph;
    std::string line;
    bool isInTable = false;
    while (std::getline(in, line)) {
        if (line.rfind("$$SOE", 0) == 0) {
            isInTable = true;
            continue;
        }
        if (line.rfind("$$EOE", 0) == 0) {
            break;
        }
        if (!isInTable) {
            continue;
        }
        // JDUT, , , longitude, latitude, delta (au), deldot,
        const std::vector<std::string> fields = splitCsv(line);
        if (fields.size() < 6) {
            throw std::runtime_error(path + " : ligne Horizons invalide « " + line + " »");
        }
        eph.t.push_back((std::stod(fields[0]) - UNIX_EPOCH_JD) * SECONDS_PER_DAY);
        eph.carrLonDeg.push_back(std::stod(fields[3]));
        eph.carrLatDeg.push_back(std::stod(fields[4]));
        eph.rAu.push_back(std::stod(fields[5]));
    }
    if (eph.t.size() < 2) {
        throw std::runtime_error(path + " : table Horizons ($$SOE ... $$EOE) absente ou trop courte");
    }
    return eph;
}

Position interpolatePosition(const Ephemeris& eph, double t) {
    if (t < eph.t.front() || t > eph.t.back()) {
        throw std::runtime_error("instant " + formatIsoTime(t) + " hors de l'éphéméride");
    }
    const size_t k = std::max<size_t>(
        1, static_cast<size_t>(std::lower_bound(eph.t.begin(), eph.t.end(), t) - eph.t.begin()));
    const double w = (t - eph.t[k - 1]) / (eph.t[k] - eph.t[k - 1]);
    const double dLon = std::remainder(eph.carrLonDeg[k] - eph.carrLonDeg[k - 1], 360.0);
    return {wrap360(eph.carrLonDeg[k - 1] + w * dLon),
            eph.carrLatDeg[k - 1] + w * (eph.carrLatDeg[k] - eph.carrLatDeg[k - 1]),
            eph.rAu[k - 1] + w * (eph.rAu[k] - eph.rAu[k - 1])};
}

std::vector<std::array<double, 3>> hourlyMeans(const VectorSeries& series, double t0,
                                               int nHours) {
    std::vector<std::array<double, 3>> sums(nHours, {0.0, 0.0, 0.0});
    std::vector<int> counts(nHours, 0);
    for (size_t i = 0; i < series.t.size(); ++i) {
        const auto& v = series.v[i];
        const long k = static_cast<long>(std::floor((series.t[i] - t0) / 3600.0));
        const bool isFinite = std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
        if (k < 0 || k >= nHours || !isFinite) {
            continue;
        }
        for (int c = 0; c < 3; ++c) {
            sums[k][c] += v[c];
        }
        ++counts[k];
    }
    std::vector<std::array<double, 3>> means(nHours, {NaN, NaN, NaN});
    for (int k = 0; k < nHours; ++k) {
        if (counts[k] > 0) {
            means[k] = {sums[k][0] / counts[k], sums[k][1] / counts[k], sums[k][2] / counts[k]};
        }
    }
    return means;
}

double parkerAlignedField(double bR, double bT, double rAu, double vR) {
    // ponytail: spirale équatoriale (sin(colatitude) = 1) ; la sonde reste à |lat| < 10°
    const double psi = std::atan(carringtonRateRadPerSecond() * rAu * AU_KM / vR);
    return bR * std::cos(psi) - bT * std::sin(psi);
}

double sourceSurfaceLongitude(double carrLonDeg, double rAu, double vR, double rssSolarRadii) {
    const double travelSeconds = (rAu * AU_KM - rssSolarRadii * SOLAR_RADIUS_KM) / vR;
    return wrap360(carrLonDeg + CARRINGTON_RATE_DEG_PER_DAY * travelSeconds / SECONDS_PER_DAY);
}

std::vector<int> sectorPolarity(const std::vector<int>& hourlySign, int windowHours) {
    const long half = windowHours / 2;
    const long n = static_cast<long>(hourlySign.size());
    std::vector<int> sector(hourlySign.size(), 0);
    for (long k = 0; k < n; ++k) {
        int sum = 0;
        for (long j = std::max(0L, k - half); j <= std::min(n - 1, k + half); ++j) {
            sum += hourlySign[j];
        }
        sector[k] = (sum > 0) - (sum < 0);
    }
    return sector;
}

std::vector<size_t> polarityChanges(const std::vector<int>& polarity) {
    std::vector<size_t> changes;
    int last = 0;
    for (size_t k = 0; k < polarity.size(); ++k) {
        if (polarity[k] == 0) {
            continue;
        }
        if (last != 0 && polarity[k] != last) {
            changes.push_back(k);
        }
        last = polarity[k];
    }
    return changes;
}
