#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include "func.h"


using namespace std;


int main(int argc, char *argv[]){
    int n, m, k, i, t, j, p = 0;
    char *name = NULL;
    double e, invar1 = 0, invar2 = 0;
    double *a = NULL;
    double *values = NULL;

    if(argc < 5 || argc > 6){
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
    if (!(convent3 >> e)){
        cout << "Error with argument e!" << endl;
        return -1;
    }
    stringstream convent4(argv[4]);
    if (!(convent4 >> k)){
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
    cout << "e - " << e << endl;
    cout << "k - " << k << endl;
    if (k == 0){
        name = argv[5];
        if (name == NULL){
            cout << "You didn't enter file name!" << endl;
            return -1;
        } else
            cout << "File name: " << name << endl;
    }
    cout << endl;

    a = (double*)malloc(n * n * sizeof(double));
    values = (double*)malloc(n * sizeof(double));

    if (a == NULL || values == NULL){
	    cout << "Not enough memory!" << endl;
	    free(a);
	    free(values);
	    return -1;
    }

    if (k == 0){
        if (fileInputMatrix(n, name, a) == -1){
            cout << "There are less number of numbers in the file!" << endl;
            free(a);
            free(values);
            return -1;
        } else if(fileInputMatrix(n, name, a) == -2){
            cout << "File did not open!" << endl;
            free(a);
            free(values);
            return -2;
        } else if(fileInputMatrix(n, name, a) == -3){
            cout << "Matrix is not symmetric!" << endl;
            free(a);
            free(values);
            return -3;
        }
    } else
        InputMatrix(k, n, a);

    cout << "Matrix A is: " << endl;
    OutputMatrix(n, n, m, a);
    cout << endl;

    for (i = 0; i < n; ++i){
        invar1 -= a[i * n + i];
        for (j = 0; j < n; ++j)
            invar2 -= a[i * n + j] * a[j * n + i];
	}

    t = clock();
    Values(n, a, values, e);
    t = clock() - t;

    for (i = 0; i < n; ++i){
        invar1 += values[i];
        invar2 += values[i] * values[i];
    }

    p = LU(n, a);

    cout << "The vector of values is:" << endl;
    OutputMatrix (1, n, m, values);
    cout << endl;

    cout << "Sum(x_i) - Sum(a_i) = " << invar1 << endl;
    cout << endl;

    cout << "Sum(x_i ^ 2) - Sum(a_ij * a_ji) = " << invar2 << endl;
    cout << endl;

    cout << "The number of sign changes of the main minors of the matrix: " << p << endl;
    cout << endl;

    cout << "Time spent:" << (double)t / CLOCKS_PER_SEC << " s" << endl;
    cout << endl;

    free(a);
    free(values);

    return 0;
}
