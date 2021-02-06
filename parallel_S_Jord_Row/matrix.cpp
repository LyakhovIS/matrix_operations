#include <iostream>
#include <fstream>
#include <cmath>
#include "func.h"


using namespace std;


void InputMatrix (int k, int n, double *a){
    int i, j;
    for(i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            a[i * n + j] = Enteringfunction(n, k, i, j);
}


double Enteringfunction(int n, int k, int i, int j){
    double element;
    if (k == 1)
        return element = n - max(i, j) + 1;
    else if (k == 2)
        return element = max(i, j);
    else if (k == 3)
        return element = abs(i - j);
    else if (k == 4)
        return element = 1.0 / (i + j + 1);
    else if (k == 5){
        if (i != (n - 1))
            return element = abs(i - j);
        else if (i == (n - 1))
            return element = 0;
    }
    return -1;
}


int fileInputMatrix(int n, char *name, double *a){
    int c = 0, i, j;
    ifstream fin(name);
    if (fin.is_open()){
        c = 0;
        for (i = 0; i < n; ++i)
            for(j = 0; j < n; ++j){
                if(!fin.eof()){
                    fin >> a[i * n + j];
                    c++;
                } else
                    break;
            }
        if (c < n*n)
            return -1;
        fin.close();
    } else
        return -2;
    return 0;
}


void InputVector (int n, double *a, double *b){
    int i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < (n + 1) / 2; ++j)
            b[i] += a[i * n + 2 * j];
}


void OutputMatrix (int l, int n, int m, double *a){
    int en, i, j;
    en = min(m, n);
    if (l == n){
        for (i = 0; i < en; ++i){
            for (j = 0; j < en; ++j)
                cout << scientific << a[i * n + j] << " ";
            cout << endl;
        }
    } else
        for (i = 0; i < en; ++i)
            cout << scientific << a[i] << endl;
}
