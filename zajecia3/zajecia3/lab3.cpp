#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double f(double x) {
    return 1.0 / (1.0 + pow(x, 3));
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

// Interpolacja Lagrangea 
double lagrange_interpolacja(const vector<double>& wezly, const vector<double>& wartosci, double x) {
    double wynik = 0.0;
    int n = wezly.size();

    for (int i = 0; i < n; ++i) {
        // Jeśli x pokrywa się dokładnie z którymś z węzłów, zwracamy jego wartość
        if (fabs(x - wezly[i]) < 1e-12) {
            return wartosci[i];
        }

        // Liczenie wielomianu Lagrange’a
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

    // Liczba węzłów (n+1 punktów)
    const int n = 10;

    // Liczba punktów do wyświetlenia na konsoli
    const int liczba_punktow = 11;

    // Generujemy węzły i obliczamy wartości funkcji f(x) w tych węzłach
    auto wezly_rown = rownoodlegle_wezly(n, a, b);
    auto wezly_czeb = czebyszew_wezly(n, a, b);

    vector<double> wartosci_rown, wartosci_czeb;
    for (auto x : wezly_rown) wartosci_rown.push_back(f(x));
    for (auto x : wezly_czeb) wartosci_czeb.push_back(f(x));

    cout << "x\tf(x)\tRównoodległe\tCzebyszew\n";
    cout << fixed << setprecision(6);

    // Wyświetlanie wyników interpolacji w równych odstępach (11 punktów)
    for (int i = 0; i < liczba_punktow; ++i) {
        double x = a + i * (b - a) / (liczba_punktow - 1); // Równe odstępy
        double y_f = f(x);
        double y_rown = lagrange_interpolacja(wezly_rown, wartosci_rown, x);
        double y_czeb = lagrange_interpolacja(wezly_czeb, wartosci_czeb, x);

        cout << x << "\t" << y_f << "\t" << y_rown << "\t";

        // Zabezpieczenie przed NaN/Inf w Czebyszewie
        if (isnormal(y_czeb) || y_czeb == 0) cout << y_czeb;
        else cout << "N/A";

        cout << endl;
    }

    ofstream f_dane("dane.txt");  // funkcja f(x)
    ofstream f_rown("rown.txt");  // interpolacja z węzłami równoodległymi
    ofstream f_czeb("czeb.txt");  // interpolacja z węzłami Czebyszewa

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
