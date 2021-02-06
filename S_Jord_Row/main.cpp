#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include "func.h"


using namespace std;


int main(int argc, char *argv[]){
    int n, m, k, i, t, j;
    int *mass = NULL;
    char *name = NULL;
    double norm = 0, maxnorm = 0;
    double *a = NULL;
    double *b = NULL;
    double *x = NULL;
    double *answer = NULL;

    cout << endl;

    if(argc < 4 || argc > 5){
        cout << "Wrong number of arguments!" << endl;
        return -1;
    }

    stringstream convent1(argv[1]);
    if (!(convent1 >> n)){
        cout << "Error with argument n!" << endl;
        return -1;
    }
    stringstream convent2(argv[2]);
    if (!(convent2 >> m)){
        cout << "Error with argument m!" << endl;
        return -1;
    }
    stringstream convent3(argv[3]);
    if (!(convent3 >> k)){
        cout << "Error with argument k!" << endl;
        return -1;
    }

    if (n < 1){
        cout << "Wrong argument n!" << endl;
        return -1;
    }

    if (m < 1){
        cout << "Wrong argument m!" << endl;
        return -1;
    }

    if (k != 0 && k != 1 && k != 2 && k != 3 && k != 4){
        cout << "Wrong argument k!" << endl;
        return -1;
    }

    cout << "Your arguments are:" << endl;
    cout << "n - " << n << endl;
    cout << "m - " << m << endl;
    cout << "k - " << k << endl;
    if (k == 0){
        name = argv[4];
        if (name == NULL){
            cout << "You didn't enter file name!" << endl;
            return -1;
        } else
            cout << "File name: " << name << endl;
    }
    cout << endl;

    a = (double*)malloc(n * n * sizeof(double));
    b = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    answer = (double*)malloc(n * sizeof(double));
    mass = (int*)malloc(n * sizeof(int));

    if (a == NULL || b == NULL || x == NULL || mass == NULL || answer == NULL){
	    cout << "Not enough memory!" << endl;
	    free(a);
        free(b);
        free(x);
        free(mass);
        free(answer);
        return -1;
    }

    for (i = 0; i < n; i++)
        mass[i] = i;

    if (k == 0){
        if (fileInputMatrix(n, name, a) == -1){
            cout << "There are less number of numbers in the file!" << endl;
            free(a);
            free(b);
            free(x);
            free(mass);
            free(answer);
            return -1;
        } else if(fileInputMatrix(n, name, a) == -2){
            cout << "File did not open!" << endl;
            free(a);
            free(b);
            free(x);
            free(mass);
            free(answer);
            return -2;
        }
    } else
        InputMatrix(k, n, a);

    InputVector(n, a, b);
    InputVector(n, a, x);

    cout << "Matrix A is: " << endl;
    OutputMatrix(n, n, m, a);
    cout << endl;

    cout << "Vector b is: " << endl;
    OutputMatrix(1, n, m, b);
    cout << endl;

    for (i = 0; i < n; ++i)
        maxnorm += abs(a[i]);

    for (i = 0; i < n; ++i){
        for (j = 0; j < n; ++j)
            norm += abs(a[i * n + j]);
        if (norm > maxnorm)
            maxnorm = norm;
        norm = 0;
    }

    if (maxnorm < 1e-10){
        for (i = 0; i < n; ++i){
            for (j = 0; j < n; ++j)
                a[i * n + j] /= maxnorm;
            x[i] /= maxnorm;
        }
    }

    t = clock();
    if (SolveSystem(n, mass, a, x, answer) == -1){
        cout << "There are no solutions: matrix A is degenerate!" << endl;
        free(a);
        free(b);
        free(x);
        free(mass);
        free(answer);
        return -1;
    }
    t = clock() - t;

    cout << "The answer is:" << endl;
    OutputMatrix (1, n, m, answer);
    cout << endl;

    if (k == 0)
        fileInputMatrix(n, name, a);
    else
        InputMatrix(k, n, a);

    cout << "||Ax - b||/||b|| = " << SolutionError(n, a, b, answer) << endl;
    cout << endl;

    cout << "The solution accuracy is: " << SolutionAccuracy(n, answer) << endl;
    cout << endl;

    cout << "Time spent:" << (double)t / CLOCKS_PER_SEC << " s" << endl;
    cout << endl;

    free(a);
    free(b);
    free(x);
    free(mass);
    free(answer);

    return 0;
}


double SolutionError(int n, double* a, double* b, double* x){
    int i, j;
    double tmp1, tmp2 = 0.0, rezult = 0.0;
    for (i = 0; i < n; ++i){
        tmp1 = 0.0;
        for (j = 0; j < n; ++j)
            tmp1 += a[i * n + j] * x[j];
        tmp1 -= b[i];
        tmp2 += b[i] * b[i];
        rezult += tmp1 * tmp1;
    }
    return sqrt(rezult / tmp2);
}


double SolutionAccuracy(int n, double* x){
    int i;
    double tmp, rezult = 0.0;
    for (i = 0; i < n; ++i){
        tmp = x[i];
        if (i % 2 == 0)
            tmp -= 1;
        rezult += tmp * tmp;
    }
    return sqrt(rezult);
}
