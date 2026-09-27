// Phase 5 : chronomètre l'ajustement étape par étape. Une ligne CSV par répétition :
//   solveur,lmax,pixels,threads,etape1,etape2,etape3,etape4,total   (secondes)
// Étapes du solveur dense    : matrice, QR + SVD, courbe en L, solution.
// Étapes du solveur par anneaux : Fourier, blocs, courbe en L, solutions.
// Usage : bench_fit <carte.fits[.gz]> <dense|rings> <lmax> <|latitude| max (°)> <répétitions>
// Threads : OMP_NUM_THREADS (noyaux Fortran) ; ceux de la BLAS se règlent à part
// (VECLIB_MAXIMUM_THREADS sur macOS, OPENBLAS_NUM_THREADS ou MKL_NUM_THREADS ailleurs).
#include "magnetosun_fortran.hpp"
#include "sh_model.hpp"
#include "synoptic_map.hpp"

#include <chrono>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>

namespace {

FitSolver parseSolver(const std::string& name) {
    if (name == "dense") {
        return FitSolver::Dense;
    }
    if (name == "rings") {
        return FitSolver::Rings;
    }
    throw std::runtime_error("solveur « " + name + " » inconnu : dense ou rings");
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr,
                     "usage : %s <carte.fits[.gz]> <dense|rings> <lmax> <|latitude| max (°)> "
                     "<répétitions>\n",
                     argv[0]);
        return 2;
    }
    try {
        const std::string solverName = argv[2];
        const FitSolver solver = parseSolver(solverName);
        const int lmax = std::stoi(argv[3]);
        const double maxAbsLatDeg = std::stod(argv[4]);
        const int repetitions = std::stoi(argv[5]);
        if (repetitions < 1) {
            throw std::runtime_error("au moins une répétition");
        }
        const SynopticMap map = readSynopticMap(argv[1]);
        std::printf("# solveur,lmax,pixels,threads,etape1,etape2,etape3,etape4,total\n");
        for (int r = 0; r < repetitions; ++r) {
            const auto start = std::chrono::steady_clock::now();
            const ShFit fit = fitSynopticMap(map, lmax, maxAbsLatDeg, -1.0, solver);
            const std::chrono::duration<double> total = std::chrono::steady_clock::now() - start;
            const auto& t = fit.stageSeconds;
            std::printf("%s,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n", solverName.c_str(), lmax,
                        fit.nPixels, ms_max_threads(), t[0], t[1], t[2], t[3], total.count());
            std::fflush(stdout);
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erreur : %s\n", e.what());
        return 1;
    }
    return 0;
}
