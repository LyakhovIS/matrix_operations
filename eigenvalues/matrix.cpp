#include <iostream>
#include <fstream>
#include <cmath>
#include "func.h"


using namespace std;
const double EPS = 1e-13;


void InputMatrix (int k, int n, double *a){
    int i, j;
    for(i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            a[i * n + j] = Enteringfunction(n, k, i, j);
}


double Enteringfunction(int n, int k, int i, int j){
    double element;
    if (k == 1){
        return element = n - max(i, j);
    }
    else if (k == 2){
        if (abs(i - j) < EPS)
            element = 2;
        else if (abs(abs(i - j) - 1) < EPS)
            element = -1;
        else
            element = 0;
        return element;
    }
    else if (k == 3){
        if (abs(i - j) < EPS && i < n && j < n)
            element = 1;
        else if (abs(j + 1 - n) < EPS)
            element = i + 1;
        else if (abs(i + 1 - n) < EPS)
            element = j + 1;
        else
            element = 0;
        if (abs(i + 1 - n) < EPS && abs(j + 1 - n) < EPS)
            element = n;
        return element;
    }
    else if (k == 4){
        return element = 1.0 / (i + j + 1);
    }
    else{
        return -1;
    }
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
    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j){
            if (abs(a[i * n + j] - a[j * n + i]) > EPS)
                return -3;
        }
    }
    return 0;
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

