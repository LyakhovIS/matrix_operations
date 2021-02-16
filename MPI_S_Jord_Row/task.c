#include <mpi.h>
#include <math.h>
#include "func.h"

const double EPS = 1e-10;

int SolveSystem(int n, int* mass, double* a, double* b, double* x, double* y, double* answer, int my_rank, int p){
    int i, j, k, l, rank, first, last, flag[1];
    int first_row, last_row, max_rows;
    double tmp, maxi, norm = 0.0, maxnorm = 0.0;
    flag[0] = 0;

    MPI_Status status;

    first_row = (n * my_rank) / p;
    last_row = (n * (my_rank + 1)) / p - 1;
    max_rows = last_row - first_row + 1;

    for (i = 0; i < n; ++i)
        maxnorm += fabs(a[i]);

    for (i = 1; i < max_rows; ++i){
        norm = 0;
        for (j = 0; j < n; ++j)
            norm += fabs(a[i * n + j]);
        if (norm > maxnorm)
            maxnorm = norm;
    }
    MPI_Allreduce (&maxnorm, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (tmp < EPS){
        for (i = 0; i < max_rows; ++i){
            for (j = 0; j < n; ++j)
                a[i * n + j] /= tmp;
            b[i] /= tmp;
        }
    }

    for (i = 0; i < n; ++i){
        if (i >= first_row && i <= last_row){
            rank = my_rank;
            maxi = fabs(a[(i - first_row) * n + i]);
            l = i;
            for (j = i + 1; j < n; ++j)
                if (fabs(a[(i - first_row) * n + j]) > maxi){
                    l = j;
                    maxi = fabs(a[(i - first_row) * n + j]);
                }

            for (k = 0; k < p; ++k){
                if (k != my_rank)
                    MPI_Send(&my_rank, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
            }
            if (fabs(maxi) < EPS){
                flag[0] = -1;
                MPI_Bcast(&flag[0], 1, MPI_INT, rank, MPI_COMM_WORLD);
                return -1;
            }
            MPI_Bcast(&flag[0], 1, MPI_INT, rank, MPI_COMM_WORLD);

            for (k = 0; k < p; ++k)
                if (k != my_rank)
                    MPI_Send(&l, 1, MPI_INT, k, 0, MPI_COMM_WORLD);

            if (l != i){
                for (j = 0; j < max_rows; ++j){
                    tmp = a[j * n + i];
                    a[j * n + i] = a[j * n + l];
                    a[j * n + l] = tmp;
                }
                k = mass[i];
                mass[i] = mass[l];
                mass[l] = k;
            }

            maxi = a[(i - first_row) * n + i];
            for (j = i, k = 0; j < n; ++j, ++k){
                a[(i - first_row) * n + j] /= maxi;
                x[k] = a[(i - first_row) * n + j];
            }
            b[i - first_row] /= maxi;
            x[k] = b[i - first_row];
            for (k = 0; k < p; ++k)
                MPI_Send(&my_rank, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
            MPI_Bcast(x, n + 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        } else{
            MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            MPI_Bcast(&flag[0], 1, MPI_INT, rank, MPI_COMM_WORLD);
            if (flag[0] == -1)
                return -1;

            MPI_Recv(&l, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

            if (l != i){
                for (j = 0; j < max_rows; ++j){
                    tmp = a[j * n + i];
                    a[j * n + i] = a[j * n + l];
                    a[j * n + l] = tmp;
                }
                k = mass[i];
                mass[i] = mass[l];
                mass[l] = k;
            }
        }

        MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (rank != my_rank)
            MPI_Bcast(x, n + 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        for (j = first_row; j <= last_row; ++j){
            if (j != i){
                b[j - first_row] -= x[n - i] * a[(j - first_row) * n + i];
                for (k = n - 1; k >= i; --k)
                    a[(j - first_row) * n + k] -= a[(j - first_row) * n + i] * x[k - i];
            }
        }
    }

    if (my_rank == 0){
        for (i = 0; i < max_rows; ++i)
            y[i] = b[i];
        if (p > 1){
            for (k = 1; k < p; ++k){
                MPI_Recv(&first, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);
                MPI_Recv(&last, 1, MPI_INT, k, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(x, n, MPI_DOUBLE, k, 2, MPI_COMM_WORLD, &status);

                for (i = first; i <= last; ++i)
                    y[i] = x[i - first];
            }
        }
    } else{
        MPI_Send(&first_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&last_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(b, n, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }

    if (my_rank == 0)
        for (i = 0; i < n; ++i)
            answer[mass[i]] = y[i];

    return 0;
}
