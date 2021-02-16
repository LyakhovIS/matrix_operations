#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "func.h"


int min(int a, int b){
    if (a >= b)
        return b;
    else
        return a;
}


int max(int a, int b){
    if (a >= b)
        return a;
    else
        return b;
}


double Enteringfunction(int n, int k, int i, int j){
    double element;
    if (k == 1)
        return element = n - max(i, j) + 1;
    else if (k == 2)
        return element = max(i, j);
    else if (k == 3){
        if (i > j)
            element = i - j;
        else
            element = j - i;
        return element;
    } else if (k == 4)
        return element = 1.0 / (i + j + 1);
    else
        return -1;
}


void InputMatrix(int n, int k, double *a, double *b, int first_row, int max_rows){
    int i, j;
    double x = 0;
    for (i = 0; i < max_rows; ++i){
        for (j = 0; j < n; ++j){
            a[i * n + j] = Enteringfunction(n, k, first_row + i, j);
            if (!(j % 2))
                x += a[i * n + j];
        }
        b[i] = x;
        x = 0;
    }
}


int fileInputMatrix(int n, char *name, double *a, double *b, double *x, int my_rank, int p){
	FILE *IN ;
	int i, j, k, err = 0;
	int first_row, last_row;
	double z = 0;

	MPI_Status status;

    first_row = (n * my_rank) / p;
    last_row = (n * (my_rank + 1)) / p - 1;

    if (my_rank == 0){
        IN = fopen(name, "r");
        if (IN == NULL){
            err = 1;
            MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
            return -1;
        }
        for (i = 0; i < n; ++i){
            for (j = 0; j < n; ++j){
                if (fscanf(IN, "%lf", &x[i * n + j]) != 1){
                    fclose(IN);
                    err = 2;
                    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    return -2;
                }
            }
        }

        fclose(IN);

        for (j = first_row; j <= last_row; ++j){
            for (k = 0; k < n; ++k){
                a[j * n + k] = x [j * n + k];
                if (!(k % 2))
                    z += x[j * n + k];
                }
            b[j - first_row] = z;
            z = 0;
        }
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (j = 1; j < p; ++j)
            MPI_Send(x, n * n, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);

    } else{
        MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (err == 1)
            return -1;
        if (err == 2)
            return -2;

        MPI_Recv(x, n * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

        for (j = first_row; j <= last_row; ++j){
            for (k = 0; k < n; ++k){
                a[(j - first_row) * n + k] = x[j * n + k];
                if (!(k % 2))
                    z += x[j * n + k];
            }
            b[j - first_row] = z;
            z = 0;
        }
    }
	return 0;
}


void OutputMatrix(int l, int n, int m, double *a, double *x, int my_rank, int p){
	int i, j, k, en, em, em1, c;
	int first_row, last_row, max_rows;

    first_row = (n * my_rank) / p;
    last_row = (n * (my_rank + 1)) / p - 1;
    max_rows = last_row - first_row + 1;

	MPI_Status status;

	en = min(m, n);
	em = min(m, max_rows);

	if (l == n){
        if (my_rank != 0){
            for (i = 0; i < max_rows; ++i){
                for (j = 0; j < en; ++j)
                    x[i * n + j] = a[i * n + j];
            }

            MPI_Send(&max_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(x, en*max_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        } else{

            for (i = 0; i < em; ++i){
                for (j = 0; j < en; ++j)
                    printf("%10.3e ", a[i * n + j]);
                printf("\n");
            }
            c = m - em;
            for (k = 1; k < p; ++k){
                MPI_Recv(&em, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(x, en*em, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
                if (c >= 0){
                    em1 = min(em, c);
                    for (i = 0; i < em1; ++i){
                        for (j = 0; j < en; ++j)
                            printf("%10.3e ", x[i * n + j]);
                        printf("\n");
                    }
                    c -= em1;
                }
            }
        }
	} else{
        if (my_rank != 0){
            for (i = 0; i < max_rows; ++i)
                x[i] = a[i];
            MPI_Send(&max_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(x, max_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        } else{
            for (i = 0; i < em; ++i)
                printf("%10.3e\n", a[i]);
            c = m - em;
            for (k = 1; k < p; ++k){
                MPI_Recv(&em, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(x, em, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
                if (c >= 0){
                    em1 = min(em, c);
                    for (i = 0; i < em1; ++i)
                        printf("%10.3e\n", x[i]);
                    c -= em1;
                }
            }
        }
	}
}


void OutputAnswer(int m, double* answer){
    int i;
    for (i = 0; i < m; ++i)
        printf("%10.3e \n", answer[i]);
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


double SolutionError(int n, double* a, double* b, double* x, double* y, double* answer, int my_rank, int p){
    int i, j, k;
    int first_row, last_row, max_rows;
    double tmp1, tmp2 = 0.0, rezult = 0.0;

    MPI_Status status;

    first_row = (n * my_rank) / p;
    last_row = (n * (my_rank + 1)) / p - 1;
    max_rows = last_row - first_row + 1;

    if (my_rank == 0){
        for (i = 0; i < max_rows; ++i){
            tmp1 = 0.0;
            for (j = 0; j < n; ++j)
                tmp1 += a[i * n + j] * answer[j];
            tmp1 -= b[i];
            tmp2 += b[i] * b[i];
            rezult += tmp1 * tmp1;
        }
        for (k = 1; k < p; ++k){
            MPI_Recv(&max_rows, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&first_row, 1, MPI_INT, k, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(x, n*max_rows, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(y, n, MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
            for (i = 0; i < max_rows; ++i){
                tmp1 = 0.0;
                for (j = 0; j < n; ++j)
                    tmp1 += x[i * n + j] * answer[j];
                tmp1 -= y[i];
                tmp2 += y[i] * y[i];
                rezult += tmp1 * tmp1;
            }
        }
    } else{
        for (i = 0; i < max_rows; ++i){
            for (j = 0; j < n; ++j)
                x[i * n + j] = a[i * n + j];
            y[i] = b[i];
        }
        MPI_Send(&max_rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&first_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(x, n*max_rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(y, n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    if(my_rank == 0) return sqrt(rezult / tmp2);
    else return 0;
}
