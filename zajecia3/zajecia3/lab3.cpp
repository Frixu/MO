#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// Funkcje do rozwiązania z zadania 1 (POPRAWIONE)
double equation1a(double x) {
    return tanh(x) + 2 * x - 2;
}

double equation1b(double x) {
    return sinh(x) + x / 4.0 - 1;
}

// Pochodne funkcji (POPRAWIONE)
double derivative1a(double x) {
    return 1.0 / (cosh(x) * cosh(x)) + 2;
}

double derivative1b(double x) {
    return cosh(x) + 1.0 / 4.0;
}

// Metoda Picarda
void picard(double (*f)(double), double x0, double tol, int max_iter) {
    cout << "Metoda Picarda:\n";
    cout << "Iteracja | x_i         | f(x_i)      | Residuum    | Blad\n";
    cout << "--------------------------------------------------------\n";

    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double x_new = x - f(x); // Picard iteration: x_{n+1} = x_n - f(x_n)
        double residual = fabs(f(x_new));
        double error = fabs(x_new - x);

        cout << setw(8) << i << " | " << setw(12) << x_new << " | " << setw(12) << f(x_new)
            << " | " << setw(12) << residual << " | " << setw(12) << error << endl;

        if (error < tol || residual < tol) {
            cout << "Zbieżność osiągnięta po " << i + 1 << " iteracjach.\n";
            return;
        }
        x = x_new;
    }
    cout << "Osiągnięto maksymalną liczbę iteracji bez zbieżności.\n";
}

// Metoda bisekcji
void bisection(double (*f)(double), double a, double b, double tol, int max_iter) {
    cout << "Metoda bisekcji:\n";
    cout << "Iteracja | a          | b          | x_i        | f(x_i)     | Residuum    | Blad\n";
    cout << "-------------------------------------------------------------------------------\n";

    if (f(a) * f(b) >= 0) {
        cout << "Nieprawidłowy przedział [a, b].\n";
        return;
    }

    double x = a;
    for (int i = 0; i < max_iter; ++i) {
        x = (a + b) / 2;
        double residual = fabs(f(x));
        double error = (b - a) / 2;

        cout << setw(8) << i << " | " << setw(12) << a << " | " << setw(12) << b << " | " << setw(12) << x
            << " | " << setw(12) << f(x) << " | " << setw(12) << residual << " | " << setw(12) << error << endl;

        if (residual < tol || error < tol) {
            cout << "Zbieżność osiągnięta po " << i + 1 << " iteracjach.\n";
            return;
        }

        if (f(a) * f(x) < 0) {
            b = x;
        }
        else {
            a = x;
        }
    }
    cout << "Osiągnięto maksymalną liczbę iteracji bez zbieżności.\n";
}

// Metoda Newtona
void newton(double (*f)(double), double (*df)(double), double x0, double tol, int max_iter) {
    cout << "Metoda Newtona:\n";
    cout << "Iteracja | x_i         | f(x_i)      | Residuum    | Blad\n";
    cout << "--------------------------------------------------------\n";

    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double fx = f(x);
        double dfx = df(x);
        if (dfx == 0) {
            cout << "Pochodna zerowa. Metoda nie może być kontynuowana.\n";
            return;
        }
        double x_new = x - fx / dfx;
        double residual = fabs(f(x_new));
        double error = fabs(x_new - x);

        cout << setw(8) << i << " | " << setw(12) << x_new << " | " << setw(12) << f(x_new)
            << " | " << setw(12) << residual << " | " << setw(12) << error << endl;

        if (error < tol || residual < tol) {
            cout << "Zbieżność osiągnięta po " << i + 1 << " iteracjach.\n";
            return;
        }
        x = x_new;
    }
    cout << "Osiągnięto maksymalną liczbę iteracji bez zbieżności.\n";
}

// Metoda siecznych
void secant(double (*f)(double), double x0, double x1, double tol, int max_iter) {
    cout << "Metoda siecznych:\n";
    cout << "Iteracja | x_i         | f(x_i)      | Residuum    | Blad\n";
    cout << "--------------------------------------------------------\n";

    double x_prev = x0;
    double x_curr = x1;
    for (int i = 0; i < max_iter; ++i) {
        double fx_prev = f(x_prev);
        double fx_curr = f(x_curr);
        if (fx_prev == fx_curr) {
            cout << "Dzielenie przez zero. Metoda nie może być kontynuowana.\n";
            return;
        }
        double x_next = x_curr - fx_curr * (x_curr - x_prev) / (fx_curr - fx_prev);
        double residual = fabs(f(x_next));
        double error = fabs(x_next - x_curr);

        cout << setw(8) << i << " | " << setw(12) << x_next << " | " << setw(12) << f(x_next)
            << " | " << setw(12) << residual << " | " << setw(12) << error << endl;

        if (error < tol || residual < tol) {
            cout << "Zbieżność osiągnięta po " << i + 1 << " iteracjach.\n";
            return;
        }
        x_prev = x_curr;
        x_curr = x_next;
    }
    cout << "Osiągnięto maksymalną liczbę iteracji bez zbieżności.\n";
}

int main() {
    double tol = 1e-6;
    int max_iter = 100;

    // Zadanie 1a: tanh(x) + 2x - 2 = 0
    cout << "Rozwiązywanie równania 1a: tanh(x) + 2x - 2 = 0\n";
    cout << "------------------------------------------------\n";
    picard(equation1a, 0.5, tol, max_iter);
    cout << "\n";
    bisection(equation1a, 0.0, 1.0, tol, max_iter);
    cout << "\n";
    newton(equation1a, derivative1a, 0.5, tol, max_iter);
    cout << "\n";
    secant(equation1a, 0.0, 1.0, tol, max_iter);
    cout << "\n";

    // Zadanie 1b: sinh(x) + x/4 - 1 = 0
    cout << "Rozwiązywanie równania 1b: sinh(x) + x/4 - 1 = 0\n";
    cout << "------------------------------------------------\n";
    picard(equation1b, 0.5, tol, max_iter);
    cout << "\n";
    bisection(equation1b, 0.0, 1.0, tol, max_iter);
    cout << "\n";
    newton(equation1b, derivative1b, 0.5, tol, max_iter);
    cout << "\n";
    secant(equation1b, 0.0, 1.0, tol, max_iter);

    return 0;
}