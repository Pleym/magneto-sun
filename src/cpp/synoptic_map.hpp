// Carte synoptique du champ radial photosphérique B_r, sur la grille canonique
// du projet : longitude de Carrington x sinus de la latitude (pixels d'aire égale).
#pragma once

#include <string>
#include <vector>

// Rayon solaire nominal (IAU 2015, résolution B3), en cm.
constexpr double SOLAR_RADIUS_CM = 6.957e10;

struct SynopticMap {
    int nLon = 0;
    int nLat = 0;
    // B_r en gauss, rangé comme le tableau Fortran br(nLon, nLat) :
    // le pixel (iLon, iLat), indices à partir de 0, est br[iLat * nLon + iLon].
    // NaN pour les pixels manquants déclarés par la carte (HMI : pôle caché).
    std::vector<double> br;
    int carringtonRotation = 0;
    std::string observationTime;  // UT de la carte
    // Centre de la colonne 0 en fraction de pixel : 0,5 (GONG, 0,5°) ou 0 (HMI, 0°).
    double lonOffsetPixels = 0.5;

    // Centre du pixel : longitude (iLon + lonOffsetPixels) * 360/nLon, en degrés.
    double longitudeDeg(int iLon) const;
    // Centre du pixel : -1 + (iLat + 1/2) * 2/nLat (iLat = 0 au pôle sud).
    double sinLatitude(int iLat) const;
};

struct FluxBalance {
    double netFluxMx;       // intégrale de B_r sur la sphère : nulle si div B = 0
    double unsignedFluxMx;  // intégrale de |B_r| sur la sphère
};

// Lit une carte synoptique GONG (mrzqs*.fits[.gz]) ou HMI (hmi.Synoptic_Mr.*.fits),
// vérifie sa grille (longitude de Carrington x sinus de latitude) et la remet dans
// l'ordre canonique : longitudes croissantes à partir de 0°, latitudes du sud au nord.
// Lève std::runtime_error si le fichier est illisible ou ne correspond pas.
SynopticMap readSynopticMap(const std::string& path);

FluxBalance computeFluxBalance(const SynopticMap& map);

// Écrit un champ sur la grille canonique en binaire brut pour gnuplot
// (float64, longitude la plus rapide, sud -> nord).
void writeRawBinary(const std::vector<double>& values, const std::string& path);
