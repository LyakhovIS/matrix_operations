#include <iostream>
#include <cmath>
#include <pthread.h>
#include "func.h"


using namespace std;
const double EPS = 1e-10;


int SolveSystem (int n, int *mass, double *a, double *x, double *answer, int thread_num, int total_threads, int *flag){
    int q, w, j, i, g, temp, last_row, first_row;
    double maxi = 0, numb = 0;

    first_row = (n * thread_num) / total_threads;
    last_row = (n * (thread_num + 1)) / total_threads - 1;

    if (thread_num == 0){
        double norm = 0, maxnorm = 0;
        for (i = 0; i < n; ++i)
                maxnorm += abs(a[i]);

        for (i = 0; i < n; ++i){
            for (j = 0; j < n; ++j)
                norm += abs(a[i * n + j]);
            if (norm > maxnorm)
                maxnorm = norm;
            norm = 0;
        }

        if (maxnorm < EPS){
            for (i = 0; i < n; ++i){
                for (j = 0; j < n; ++j)
                    a[i * n + j] /= maxnorm;
                x[i] /= maxnorm;
            }
        }
    }
    synchronize(total_threads);

    for (q = 0; q < n; ++q){
        synchronize(total_threads);
        if (thread_num == 0){
            maxi = abs(a[q * n + q]);
            i = q;
            for (w = q + 1; w < n; ++w)
                if(abs(a[q * n + w]) > maxi){
                    i = w;
                    maxi = abs(a[q * n + w]);
                }
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

            if (abs(maxi) < EPS){
                flag[0] = -1;
            }

            maxi = a[q * n + q];
            a[q * n + q] = 1;
            for (j = q + 1; j < n; ++j)
                a[q * n + j] /= maxi;
            x[q] /= maxi;
        }
        synchronize(total_threads);
        if (flag[0] == -1)
                return -1;

        for (i = first_row; i < last_row + 1; ++i){
            if (i != q){
                x[i] -= x[q] * a[i * n + q];
                for (g = n - 1; g > -1; --g)
                    a[i * n + g] -= a[q * n + g] * a[i * n + q];
            }
        }
    }

    if (thread_num == 0){
        for (i = 0; i < n; ++i)
            answer[mass[i]] = x[i];
    }

    return 0;
}
