#include <iostream>
#include <math.h>
#include "func.h"


using namespace std;


void diag3(int n, double* a){
	int k, m, l;
	double x, y;
	for (k = 0; k < n - 2; ++k){
		x = 0;
		for (m = k + 2; m < n; ++m)
			x += a[m * n + k] * a[m * n + k];
		y = sqrt(x + a[(k + 1) * n + k] * a[(k + 1) * n + k]);
		if (y < 1e-16){
            a[(k + 1) * n + k] = 0;
            a[(k + 2) * n + k] = 0;
            continue;
        }
        if (x < 1e-16){
            a[(k + 2) * n + k] = 0;
            continue;
        }
		a[(k + 1) * n + k] -= y;
		x = 1.0 / sqrt(x + a[(k + 1) * n + k] * a[(k + 1) * n + k]);
		for (m = k + 1; m < n; ++m)
			a[m * n + k] *= x;
		for (m = k + 1; m < n; ++m){
			x = 0;
			for (l = k + 1; l < n; l++)
				x += a[m * n + l] * a[l * n + k];
			a[k * n + m] = x;
		}
		x = 0;
		for (m = k + 1; m < n; ++m)
			x += a[k * n + m] * a[m * n + k];
		x *= 2;
		for (m = k + 1; m < n; ++m)
			a[k * n + m] = 2.0 * a[k * n + m] - x * a[m * n + k];
		for (m = k + 1; m < n; ++m)
			for (l = k + 1; l < n; ++l)
				a[m * n + l] -= a[k * n + m] * a[l * n + k] + a[k * n + l] * a[m * n + k];
		a[(k + 1) * n + k] = y;
		a[k * n + (k + 1)] = y;
		for (m = k + 2; m < n; ++m){
			a[m * n + k] = 0;
			a[k * n + m] = 0;
		}
	}
}


double Norm(int n, double* a){
	int i, j;
	double tmp, norm = 0;
	for (i = 0; i < n; ++i){
		tmp = 0;
		for (j = 0; j < n; ++j)
			tmp += fabs(a[i * n + j]);
		if (norm < tmp)
			norm = tmp;
	}
	return norm;
}


int n_(int n, double* a, double lambda, double e){
	int i, number;
	double l;
	l = a[0] - lambda;

	if (l < 0){
        number = 1;
    } else{
        number = 0;
    }

	for (i = 1; i < n; ++i){
		if (fabs(l) < e)
			l = e;
		l = (a[i * n + i] - lambda) - a[i * n + (i - 1)] * a[(i - 1) * n + i] / l;
		if (l < 0)
			++number;
	}
	return number;
}


void Values(int n, double* a, double* values, double e){
	int i = 0, j, tmp;
	double x, left, right, middle, numb;
	x = Norm(n, a) + e;
	right = x;
	left = -x;
	diag3(n, a);
	while (i < n){
		while (right - left > e && (n_(n, a, right, e) - n_(n, a, left, e)) != 0){
			middle = (left + right) / 2;
            numb = right - left;
			if (n_(n, a, middle, e) < i + 1)
				left = middle;
			else
				right = middle;

            if(abs(right - left - numb) < 1e-81)
                break;
            //cout << right - left - numb << endl;
        }

        middle = (left + right) / 2;
        tmp = n_(n, a, right, e) - n_(n, a, left, e);
        //cout << tmp << endl;

        for (j = 0; j < tmp; ++j)
            values[i + j] = middle;

        i += tmp;

        right = x;
        left = middle;
    }
}

int LU (int n, double* a){
    int i, j = 0;
    double temp1, temp2;

    temp1 = a[0];
    diag3(n, a);

    for (i = 1; i < n; ++i){
        temp2 = a[i * n + i] - a[i * n + (i - 1)] * a[(i - 1) * n + i] / temp1;
        if (temp1 * temp2 < 0)
            j++;
        temp1 = temp2;
    }
    return j;
}
