// Modèle de B_r en harmoniques sphériques et champ PFSS, côté C++ : appelle les
// noyaux Fortran et lit/écrit les coefficients.
#pragma once

#include "synoptic_map.hpp"

#include <array>
#include <string>
#include <vector>

// Rayon de la surface source usuel (Altschuler & Newkirk 1969), en rayons solaires.
constexpr double DEFAULT_SOURCE_SURFACE_RADIUS = 2.5;

struct ShCoefficients {
    int lmax = 0;
    // g(l, m) en g[m * (lmax + 1) + l], comme le tableau Fortran g(0:lmax, 0:lmax) ; idem h.
    std::vector<double> g;
    std::vector<double> h;
};

struct ShFit {
    ShCoefficients coeffs;
    std::vector<double> lambdas;  // courbe en L : lambda, |A c - b|, |c|
    std::vector<double> residualNorms;
    std::vector<double> solutionNorms;
    int cornerIndex;
    double lambda;           // lambda utilisé
    double conditionNumber;  // sigma_max / sigma_min
    int nPixels;             // pixels utilisés
    // Temps des 4 étapes (s) : dense = matrice, QR + SVD, courbe en L, solution ;
    // anneaux = Fourier, blocs, courbe en L, solutions.
    std::array<double, 4> stageSeconds;
};

// Rings : anneaux de latitude complets, découpage exact par blocs (m ; cos/sin), rapide.
// Dense : pixels quelconques, grande matrice (référence, et cas hors grille régulière).
// Les deux donnent la même solution sur une carte synoptique (tests/test_sh_rings.f90).
enum class FitSolver { Rings, Dense };

// Ajuste B_r sur les pixels de |latitude| <= maxAbsLatDeg, pixels manquants exclus
// (anneaux : rangées incomplètes exclues). lambda < 0 : lambda au coin de la courbe
// en L ; sinon lambda imposé (0 : moindres carrés ordinaires).
ShFit fitSynopticMap(const SynopticMap& map, int lmax, double maxAbsLatDeg, double lambda = -1.0,
                     FitSolver solver = FitSolver::Rings);

// Coefficients de B_r au rayon r (rayons solaires), monopôle retiré.
ShCoefficients pfssCoefficients(const ShCoefficients& coeffs, double rss, double r);

// B_r du modèle PFSS au rayon r aux points (cos(colatitude), longitude en radians).
// r = 1 : B_r photosphérique du modèle, sans le monopôle.
std::vector<double> pfssBr(const ShCoefficients& coeffs, double rss, double r,
                           const std::vector<double>& cosTheta, const std::vector<double>& phi);

// Même chose aux centres des pixels d'une grille canonique nLon x nLat
// (rangement de SynopticMap::br).
std::vector<double> pfssBrOnGrid(const ShCoefficients& coeffs, double rss, double r, int nLon,
                                 int nLat);

// Fichier texte « l m g h », précédé de lignes de commentaire dont « # lmax N ».
void writeShCoefficients(const ShCoefficients& coeffs, const std::string& path,
                         const std::string& comment);
ShCoefficients readShCoefficients(const std::string& path);
