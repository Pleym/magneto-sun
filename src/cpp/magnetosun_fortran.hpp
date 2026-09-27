// Prototypes des noyaux Fortran (src/fortran/capi.f90).
// Coefficients g, h : tableaux Fortran g(0:lmax, 0:lmax), soit g(l, m) en
// g[m * (lmax + 1) + l]. Points : colatitude par cos(theta), longitude phi en radians.
// timings : temps des 4 étapes de l'ajustement, en secondes.
#pragma once

extern "C" {

// Solveur dense, pixels quelconques (moindres carrés + Tikhonov). lambdaIn < 0 : lambda
// au coin de la courbe en L (nLambda >= 3 points). info : 0 si succès, -1 si
// nPix < (lmax+1)^2, > 0 LAPACK. Étapes : matrice, QR + SVD, courbe en L, solution.
void ms_fit(int lmax, int nPix, const double* cosTheta, const double* phi, const double* br,
            int nLambda, double lambdaIn, double* g, double* h, double* lambdas,
            double* residualNorms, double* solutionNorms, int* cornerIndex, double* lambdaUsed,
            double* conditionNumber, int* info, double* timings);

// Solveur par anneaux complets de longitude régulière, même résultat que ms_fit :
// br[j * nLon + i] = pixel i de l'anneau j, cosTheta[j] par anneau, phi0 = longitude du
// premier pixel. info : -1 si nRings <= lmax, -2 si 2 lmax >= nLon, > 0 LAPACK.
// Étapes : Fourier, blocs (m ; cos/sin), courbe en L, solutions.
void ms_fit_rings(int lmax, int nLon, int nRings, double phi0, const double* cosTheta,
                  const double* br, int nLambda, double lambdaIn, double* g, double* h,
                  double* lambdas, double* residualNorms, double* solutionNorms,
                  int* cornerIndex, double* lambdaUsed, double* conditionNumber, int* info,
                  double* timings);

// Coefficients de B_r du modèle PFSS au rayon r (rayons solaires), monopôle retiré.
void ms_pfss_coefs(int lmax, const double* g, const double* h, double rss, double r, double* gr,
                   double* hr);

// B_r du modèle PFSS au rayon r aux n points.
void ms_pfss_br(int lmax, const double* g, const double* h, double rss, double r, int n,
                const double* cosTheta, const double* phi, double* values);

// Nombre de threads OpenMP des noyaux Fortran.
int ms_max_threads();
}
