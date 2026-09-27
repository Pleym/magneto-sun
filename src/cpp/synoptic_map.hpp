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
    std::vector<double> br;
    int carringtonRotation = 0;
    std::string observationTime;  // UT du dernier magnétogramme intégré

    // Centre du pixel : longitude (iLon + 1/2) * 360/nLon, en degrés.
    double longitudeDeg(int iLon) const;
    // Centre du pixel : -1 + (iLat + 1/2) * 2/nLat (iLat = 0 au pôle sud).
    double sinLatitude(int iLat) const;
};

struct FluxBalance {
    double netFluxMx;       // intégrale de B_r sur la sphère : nulle si div B = 0
    double unsignedFluxMx;  // intégrale de |B_r| sur la sphère
};

// Lit une carte synoptique horaire GONG (mrzqs*.fits ou .fits.gz), vérifie sa
// grille et la remet sur la grille canonique (première colonne à 0,5 pixel de 0°).
// Lève std::runtime_error si le fichier est illisible ou ne correspond pas.
SynopticMap readGongMap(const std::string& path);

FluxBalance computeFluxBalance(const SynopticMap& map);
