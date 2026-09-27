// Vérifie la chaîne de compilation mixte Fortran + C++ et fige la convention
// de stockage : un tableau Fortran A(m, n) est rangé par colonnes, donc
// A(i, j) (indices 1..m, 1..n) se lit a[(j-1)*m + (i-1)] côté C++.
#include <cstdio>
#include <vector>

extern "C" void fill_matrix(int m, int n, double* a);

int main() {
    const int m = 3, n = 4;
    std::vector<double> a(m * n);
    fill_matrix(m, n, a.data());

    for (int j = 1; j <= n; ++j) {
        for (int i = 1; i <= m; ++i) {
            const double got = a[(j - 1) * m + (i - 1)];
            const double expected = 10.0 * i + j;
            if (got != expected) {
                std::printf("A(%d,%d) = %g, attendu %g\n", i, j, got, expected);
                return 1;
            }
        }
    }
    std::puts("interop Fortran/C++ OK");
    return 0;
}
