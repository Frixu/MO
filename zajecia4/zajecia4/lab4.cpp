#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm> // dla funkcji max()

using namespace std;

// Funkcje układu równań
double f1(double x, double y, double z) {
    return x * x + y * y + z * z - 4.0;
}

double f2(double x, double y) {
    return x * x + (y * y) / 2 - 1.0;
}

double f3(double x, double y) {
    return x * y - 0.5;
}

// Macierz Jakobiego
void computeJacobian(double x, double y, double z, double J[3][3]) {
    J[0][0] = 2.0 * x; J[0][1] = 2.0 * y; J[0][2] = 2.0 * z;
    J[1][0] = 2.0 * x; J[1][1] = y; J[1][2] = 0.0;
    J[2][0] = y; J[2][1] = x; J[2][2] = 0.0;
}

// Rozwiązanie układu równań liniowych metodą Gaussa
void solveLinearSystem(double A[3][3], double b[3], double dx[3]) {
    // Prosta implementacja eliminacji Gaussa
    for (int i = 0; i < 3; ++i) {
        // Pivot
        double max = abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < 3; ++k) {
            if (abs(A[k][i]) > max) {
                max = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Zamiana wierszy
        for (int k = i; k < 3; ++k) {
            swap(A[maxRow][k], A[i][k]);
        }
        swap(b[maxRow], b[i]);

        // Eliminacja
        for (int k = i + 1; k < 3; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 3; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Podstawienie wsteczne
    for (int i = 2; i >= 0; --i) {
        dx[i] = b[i];
        for (int j = i + 1; j < 3; ++j) {
            dx[i] -= A[i][j] * dx[j];
        }
        dx[i] /= A[i][i];
    }
}

// Uogólniona metoda Newtona
void generalizedNewtonMethod(double x0, double y0, double z0, double tol1, double tol2, double tol3, int maxIter) {
    double x = x0, y = y0, z = z0;
    double J[3][3];
    double F[3], dx[3];
    double prev_x = x0, prev_y = y0, prev_z = z0;

    cout << "Iteracja\tx\t\ty\t\tz\t\tEstymator błędu\tResiduum" << endl;
    cout << fixed << setprecision(10);

    for (int iter = 0; iter < maxIter; ++iter) {
        // Oblicz wartości funkcji
        F[0] = -f1(x, y, z);
        F[1] = -f2(x, y);
        F[2] = -f3(x, y);

        // Oblicz macierz Jakobiego
        computeJacobian(x, y, z, J);

        // Rozwiąż układ równań liniowych
        solveLinearSystem(J, F, dx);

        // Aktualizacja rozwiązania
        x += dx[0];
        y += dx[1];
        z += dx[2];

        // Oblicz estymator błędu (maksymalna zmiana składowych)
        double error_est = max({ abs(dx[0]), abs(dx[1]), abs(dx[2]) });

        // Oblicz residuum
        double residual = sqrt(f1(x, y, z) * f1(x, y, z) + f2(x, y) * f2(x, y) + f3(x, y) * f3(x, y));

        // Wyświetl wyniki iteracji
        cout << iter + 1 << "\t" << x << "\t" << y << "\t" << z << "\t" << error_est << "\t" << residual << endl;

        // Sprawdź warunki stopu
        if (error_est < tol1) {
            cout << "Warunek stopu 1 (estymator błędu < " << tol1 << ") spełniony." << endl;
            break;
        }
        if (residual < tol2) {
            cout << "Warunek stopu 2 (residuum < " << tol2 << ") spełniony." << endl;
            break;
        }
        if (iter > 0 && abs(x - prev_x) < tol3 && abs(y - prev_y) < tol3 && abs(z - prev_z) < tol3) {
            cout << "Warunek stopu 3 (zmiana rozwiązania < " << tol3 << ") spełniony." << endl;
            break;
        }

        prev_x = x;
        prev_y = y;
        prev_z = z;
    }

    cout << "\nOstateczne rozwiązanie:" << endl;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;
    cout << "z = " << z << endl;
   /* cout << "Sprawdzenie:" << endl;
    cout << "f1(x,y,z) = " << f1(x, y, z) << endl;
    cout << "f2(x,y) = " << f2(x, y) << endl;
    cout << "f3(x,y) = " << f3(x, y) << endl;
    */
    }

int main() {
    // Początkowe przybliżenie (należy wybrać odpowiednie)
    double x0 = 0.5, y0 = 1.0, z0 = 1.0;

    // Tolerancje dla trzech kryteriów stopu
    double tol1 = 1e-8;  // estymator błędu
    double tol2 = 1e-8;  // residuum
    double tol3 = 1e-8;  // zmiana rozwiązania

    // Maksymalna liczba iteracji
    int maxIter = 50;

    generalizedNewtonMethod(x0, y0, z0, tol1, tol2, tol3, maxIter);

    return 0;
}