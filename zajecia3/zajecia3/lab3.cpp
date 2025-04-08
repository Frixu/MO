#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

void printValue(double val) {
    if (val == floor(val)) {
        cout << setw(10) << static_cast<int>(val);
    }
    else {
        cout << setw(10) << val;
    }
}

void printMatrix(const vector<vector<double>>& matrix, const string& name) {
    cout << "--- " << name << " ---" << endl;
    for (const auto& row : matrix) {
        for (double val : row) {
            printValue(val);
        }
        cout << endl;
    }
}

void printVector(const vector<double>& vec, const string& name) {
    cout << "--- " << name << " ---" << endl;
    for (double val : vec) {
        printValue(val);
    }
    cout << endl;
}

void LUDecompositionWithPartialPivoting(vector<vector<double>>& A, vector<int>& pivot) {
    int n = A.size();
    pivot.resize(n);

    for (int i = 0; i < n; ++i) {//wektor permutacji
        pivot[i] = i;
    }

    for (int k = 0; k < n; ++k) {
        int max_row = k;
        double max_val = abs(A[k][k]);

        for (int i = k + 1; i < n; ++i) {//znajdujemy maksymalny element
            if (abs(A[i][k]) > max_val) {
                max_val = abs(A[i][k]);
                max_row = i;
            }
        }

        if (max_row != k) {//zamieniamy wiersze jesli trzeba
            swap(A[k], A[max_row]);
            swap(pivot[k], pivot[max_row]);
        }

        for (int i = k + 1; i < n; ++i) { //eliminacja gaussa
            A[i][k] /= A[k][k]; //obliczamy mnoznik
            for (int j = k + 1; j < n; ++j) {
                A[i][j] -= A[i][k] * A[k][j];//aktualizujemy podmacierz
            }
        }
    }
}

vector<double> SolveWithLU(const vector<vector<double>>& LU, const vector<int>& pivot, const vector<double>& b) {
    int n = LU.size();
    vector<double> x(n), y(n);

    vector<double> pb(n);
    for (int i = 0; i < n; ++i) {
        pb[i] = b[pivot[i]];//permutacja wektora b
    }

    for (int i = 0; i < n; ++i) {
        y[i] = pb[i];
        for (int j = 0; j < i; ++j) {//Ly = pb
            y[i] -= LU[i][j] * y[j];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {//Ux = y
            x[i] -= LU[i][j] * x[j];
        }
        x[i] /= LU[i][i];
    }

    return x;
}

void extractLandU(const vector<vector<double>>& LU, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = LU.size();
    L.assign(n, vector<double>(n, 0.0));
    U.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) { //rozdziela macierz na L i U
        L[i][i] = 1.0;
        for (int j = 0; j < i; ++j) {
            L[i][j] = LU[i][j];
        }
        for (int j = i; j < n; ++j) {
            U[i][j] = LU[i][j];
        }
    }
}

int main() {
    cout << fixed;

    // Dane wejściowe z przykładu
    vector<vector<double>> A = {
        {5, 4, 3, 2, 1},
        {10, 8, 7, 6, 5},
        {-1, 2, -3, 4, -5},
        {6, 5, -4, 3, -2},
        {1, 2, 3, 4, 5}
    };

    vector<double> b = { 37, 99, -9, 12, 53 };

    vector<vector<double>> A_copy = A;
    vector<double> b_copy = b;

        cout << "ROZWIAZANIE ===" << endl;
        printMatrix(A, "MATRIX A");
        printVector(b, "VECTOR B");

        vector<int> pivot;
        LUDecompositionWithPartialPivoting(A, pivot);

        vector<vector<double>> L, U;
        extractLandU(A, L, U);

        printMatrix(L, "MATRIX L");
        printMatrix(U, "MATRIX U");

        cout << "--- PERMUTATION VECTOR ---" << endl;
        for (int val : pivot) {
            cout << setw(10) << val;
        }
        cout << endl;

        vector<double> y = SolveWithLU(A, pivot, b_copy);
        printVector(y, "VECTOR Y");

        vector<double> x = SolveWithLU(A, pivot, b_copy);
        printVector(x, "VECTOR X");

    return 0;
}
