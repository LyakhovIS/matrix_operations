#include <iostream>
#include <cmath>
#include "func.h"


using namespace std;
const double EPS = 1e-13;


int SolveSystem (int n, int *mass, double *a, double *x, double *answer){
    int q, w, j, i, g, temp;
    double maxi, numb;

    for (q = 0; q < n; ++q){
        maxi = abs(a[q * n + q]);
        i = q;
        for(w = q + 1; w < n; ++w)
            if(abs(a[q * n + w]) > maxi){
                i = w;
                maxi = abs(a[q * n + w]);
            }
        if (abs(maxi) < EPS)
            return -1;
        if (i != q){
            for (j = 0; j < n; ++j){
                numb = a[j * n + q];
                a[j * n + q] = a[j * n + i];
                a[j * n + i] = numb;
            }
            temp = mass[i];
            mass[i] = mass[q];
            mass[q] = temp;
        }
        maxi = a[q * n + q];
        a[q * n + q] = 1;
        for (j = q + 1; j < n; ++j)
            a[q * n + j] /= maxi;
        x[q] /= maxi;
        for (i = 0; i < n; ++i){
            if (i != q){
                x[i] -= x[q] * a[i * n + q];
                for (g = n - 1; g >= q; --g)
                    a[i * n + g] -= a[q * n + g] * a[i * n + q];
            }
        }
    }
    for (i = 0; i < n; ++i)
        answer[mass[i]] = x[i];
    return 0;
}
