// Phase 4 : lecture des données in situ (AMDA, JPL Horizons) et physique de la
// polarité : spirale de Parker, projection balistique, changements de secteur.
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

constexpr double AU_KM = 149597870.7;
constexpr double SOLAR_RADIUS_KM = 695700.0;  // comme SOLAR_RADIUS_CM (IAU 2015)
// Rotation du repère de Carrington (période sidérale 25,38 jours), en degrés par jour.
constexpr double CARRINGTON_RATE_DEG_PER_DAY = 360.0 / 25.38;

// Temps en secondes depuis 1970-01-01T00:00:00 UTC.
double parseIsoTime(const std::string& text);  // « AAAA-MM-JJTHH:MM:SS[.sss] »
std::string formatIsoTime(double t);           // « AAAA-MM-JJTHH:MM »

struct VectorSeries {
    std::vector<double> t;
    std::vector<std::array<double, 3>> v;
};

// Fichier ASCII d'AMDA : lignes « date x y z » (NaN acceptés), commentaires en « # ».
VectorSeries readAmdaVectors(const std::string& path);

struct Ephemeris {
    std::vector<double> t;
    std::vector<double> carrLonDeg;
    std::vector<double> carrLatDeg;
    std::vector<double> rAu;
};

// Table d'observation JPL Horizons en CSV (cible : le Soleil, observateur : la sonde,
// quantités 14 et 20, dates en JD) : le point sous la sonde donne sa longitude et sa
// latitude de Carrington, « delta » sa distance au Soleil.
Ephemeris readHorizons(const std::string& path);

struct Position {
    double carrLonDeg;
    double carrLatDeg;
    double rAu;
};

// Interpolation linéaire (longitude déroulée) ; lève une erreur hors de l'éphéméride.
Position interpolatePosition(const Ephemeris& ephemeris, double t);

// Moyennes des échantillons finis sur [t0 + k h, t0 + (k + 1) h) ; NaN si aucun.
std::vector<std::array<double, 3>> hourlyMeans(const VectorSeries& series, double t0, int nHours);

// Composante de B le long de la spirale de Parker sortante nominale :
// B_R cos(psi) - B_T sin(psi), avec tan(psi) = Omega r / v_R. Positive : polarité sortante.
double parkerAlignedField(double bR, double bT, double rAu, double vR);

// Longitude de Carrington (°) d'où est parti, sur la surface source, le vent mesuré en r :
// il a voyagé tau = (r - Rss) / v pendant que le Soleil tournait de Omega tau.
double sourceSurfaceLongitude(double carrLonDeg, double rAu, double vR, double rssSolarRadii);

// Polarité de secteur : signe de la somme des polarités horaires (+1, -1, 0 inconnue)
// sur une fenêtre centrée de windowHours heures ; efface les inversions brèves.
std::vector<int> sectorPolarity(const std::vector<int>& hourlySign, int windowHours);

// Heures où la polarité change par rapport à la dernière heure de polarité connue.
std::vector<size_t> polarityChanges(const std::vector<int>& polarity);
