#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

// Definicja π jeśli nie jest zdefiniowana w systemie
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Poprawiona funkcja Rungego: f(x) = 1 / (1 + 25x^2)
double f(double x) {
    return 1.0 / (1.0 + 25 * x * x);
}

// Generowanie węzłów równoodległych na przedziale [a, b]
vector<double> rownoodlegle_wezly(int n, double a, double b) {
    vector<double> wezly(n + 1);
    double h = (b - a) / n;
    for (int i = 0; i <= n; ++i) {
        wezly[i] = a + i * h;
    }
    return wezly;
}

// Generowanie węzłów Czebyszewa na przedziale [a, b]
vector<double> czebyszew_wezly(int n, double a, double b) {
    vector<double> wezly(n + 1);
    for (int i = 0; i <= n; ++i) {
        double theta = (2 * i + 1) * M_PI / (2 * (n + 1));
        wezly[i] = (a + b) / 2 + (b - a) / 2 * cos(theta); // skalowanie na [a, b]
    }
    return wezly;
}

// Interpolacja Lagrange’a
double lagrange_interpolacja(const vector<double>& wezly, const vector<double>& wartosci, double x) {
    double wynik = 0.0;
    int n = wezly.size();

    for (int i = 0; i < n; ++i) {
        if (fabs(x - wezly[i]) < 1e-12) {
            return wartosci[i]; // Unikamy dzielenia przez 0 jeśli x pokrywa się z węzłem
        }

        double li = 1.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                li *= (x - wezly[j]) / (wezly[i] - wezly[j]);
            }
        }
        wynik += wartosci[i] * li;
    }

    return wynik;
}

int main() {
    // Zakres interpolacji
    const double a = -1.0;
    const double b = 1.0;

    // Liczba węzłów: n+1
    const int n = 10;

    // Liczba punktów do wypisania na konsolę
    const int liczba_punktow = 11;

    // Węzły i wartości funkcji w tych węzłach
    auto wezly_rown = rownoodlegle_wezly(n, a, b);
    auto wezly_czeb = czebyszew_wezly(n, a, b);

    vector<double> wartosci_rown, wartosci_czeb;
    for (auto x : wezly_rown) wartosci_rown.push_back(f(x));
    for (auto x : wezly_czeb) wartosci_czeb.push_back(f(x));

    // Wypisanie nagłówka
    cout << "x\tf(x)\tRównoodległe\tCzebyszew\n";
    cout << fixed << setprecision(6);

    // Obliczenie i wypisanie wartości interpolowanych w równych odstępach
    for (int i = 0; i < liczba_punktow; ++i) {
        double x = a + i * (b - a) / (liczba_punktow - 1);
        double y_f = f(x);
        double y_rown = lagrange_interpolacja(wezly_rown, wartosci_rown, x);
        double y_czeb = lagrange_interpolacja(wezly_czeb, wartosci_czeb, x);

        cout << x << "\t" << y_f << "\t" << y_rown << "\t";

        if (isnormal(y_czeb) || y_czeb == 0) cout << y_czeb;
        else cout << "N/A";

        cout << endl;
    }

    // Tworzenie plików do wykresu (100 punktów)
    ofstream f_dane("dane.txt");
    ofstream f_rown("rown.txt");
    ofstream f_czeb("czeb.txt");

    for (int i = 0; i < liczba_punktow * 10; ++i) {
        double x = a + i * (b - a) / (liczba_punktow * 10 - 1);
        f_dane << x << " " << f(x) << endl;
        f_rown << x << " " << lagrange_interpolacja(wezly_rown, wartosci_rown, x) << endl;
        f_czeb << x << " " << lagrange_interpolacja(wezly_czeb, wartosci_czeb, x) << endl;
    }

    f_dane.close();
    f_rown.close();
    f_czeb.close();

    return 0;
}
