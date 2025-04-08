#include <iostream>
#include <vector>

using namespace std;

// Struktura do przechowywania przetworzonych elementów macierzy A
struct ThomasLU {
    vector<double> gamma; // przetworzona główna przekątna
    vector<double> u;     // przetworzona nadprzekątna
    vector<double> l;     // dolna przekątna (niezmieniona)
};

// Procedura 1: działa tylko na macierzy A
ThomasLU processMatrixA(const vector<double>& d, const vector<double>& l, const vector<double>& u) {
    int n = d.size();
    vector<double> gamma(n);
    vector<double> u_mod = u;

    gamma[0] = d[0];//pierwszego elementu nie zmieniamy

    for (int i = 1; i < n; ++i) {
        u_mod[i - 1] = u[i - 1]; // górna przekątna robiona zeby nie zmieniac orginalu
        gamma[i] = d[i] - l[i - 1] * u_mod[i - 1] / gamma[i - 1];//gamma-zmodifikowana górna przekątna
    }

    return ThomasLU{ gamma, u_mod, l };
}

// Procedura 2: działa tylko na wektorze b
vector<double> processVectorB(const ThomasLU& A, const vector<double>& b) {
    int n = b.size();
    vector<double> r(n); //wektor pomocniczy do przechowywania r
    vector<double> x(n);

    // Eliminacja
    r[0] = b[0];//pierwszy element zostaje taki sam

    for (int i = 1; i < n; ++i) {
        r[i] = b[i] - A.l[i - 1] * r[i - 1] / A.gamma[i - 1];
    }

    // Podstawienie wsteczne
    x[n - 1] = r[n - 1] / A.gamma[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = (r[i] - A.u[i] * x[i + 1]) / A.gamma[i];
    }

    return x;
}

int main() {
    // Dane z zadania
    vector<double> d = { 100, 200, 300, 200, 100 };     // przekątna główna
    vector<double> l = { 2, 4, -6, -8 };              // dolna przekątna
    vector<double> u = { -1, -3, 5, -7 };              // nadprzekątna
    vector<double> b = { 199, 195, 929, 954, 360 };    // wektor b

    // Procedura 1 – tylko A
    ThomasLU A_transformed = processMatrixA(d, l, u);

    // Procedura 2 – tylko b
    vector<double> x = processVectorB(A_transformed, b);

    // Wynik
    cout << "Rozwiązanie x:\n";
    for (double xi : x)
        cout << xi << endl;

    return 0;
}
