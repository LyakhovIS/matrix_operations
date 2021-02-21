#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "func.h"


int main(int argc, char *argv[]){
    int n, k, m, p, my_rank, i, err1 = 0, err2;
    int first_row, last_row, max_rows;
    double t, error = 1.0;
    char *name;
    int *mass = NULL;
    double *a = NULL;
    double *b = NULL;
    double *x = NULL;
    double *y = NULL;
    double *answer = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) printf("\n");

	if ((argc != 5) && (argc != 4)) {
        if (my_rank == 0) printf("Wrong number of arguments!\n");
        MPI_Finalize();
        return -1;
	}

    n = atoi(argv[1]);
    m = atoi(argv[2]);
    k = atoi(argv[3]);

    if (n < 1){
        if (my_rank == 0) printf("Wrong argument n!\n");
        MPI_Finalize();
        return -2;
    }

    if (m < 1 || m > n){
        if (my_rank == 0) printf("Wrong argument m!\n");
        MPI_Finalize();
        return -2;
    }

    if (k != 0 && k != 1 && k != 2 && k != 3 && k != 4){
        if (my_rank == 0) printf("Wrong argument k!\n");
        MPI_Finalize();
        return -2;
    }

    if(p > n || p < 1){
        if (my_rank == 0) printf("Wrong number of threads!\n");
        MPI_Finalize();
        return -2;
    }


    if (my_rank == 0) printf("Your arguments are:\n");
    if (my_rank == 0) printf("n - %d\n", n);
    if (my_rank == 0) printf("m - %d\n", m);
    if (my_rank == 0) printf("k - %d\n", k);
    if (my_rank == 0) printf("p - %d\n", p);
    if (my_rank == 0){
        if (k == 0){
            name = argv[4];
            if (name == NULL)
                err1 = 1;
             else
                if (my_rank == 0) printf("File name: %s\n", name);
        }
    }

    if (my_rank == 0) printf("\n");

    MPI_Bcast(&err1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (err1 == 1) {
		if (my_rank == 0) printf("\nYou didn't enter file name!\n");
		MPI_Finalize();
		return -3;
	}


    first_row = (n * my_rank) / p;
    last_row = (n * (my_rank + 1)) / p - 1;
    max_rows = last_row - first_row + 1;

    a = (double*)malloc(n * max_rows * sizeof(double));
    b = (double*)malloc(n * max_rows * sizeof(double));
    x = (double*)malloc(n * n * sizeof(double));
    y = (double*)malloc(n * sizeof(double));
    answer = (double*)malloc(n * sizeof(double));
    mass = (int*)malloc(n * sizeof(int));

    if (a == NULL || b == NULL || x == NULL || y == NULL || answer == NULL || mass == NULL){
        err1 = 1;
    }

    MPI_Allreduce(&err1, &err2, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (err2 == 1) {
        if (my_rank == 0) printf("\nNot enough memory!\n");
        free(a);
        free(b);
        free(x);
        free(y);
        free(answer);
        free(mass);
        MPI_Finalize();
        return -4;
    }

    if (k == 0){
        i = fileInputMatrix(n, argv[4], a, b, x, my_rank, p);
        if (i == -1){
            if (my_rank == 0) printf("\nFile did not open!\n");
            free(a);
            free(b);
            free(x);
            free(y);
            free(mass);
            free(answer);
            MPI_Finalize();
            return -5;
         } else if(i == -2){
            if (my_rank == 0) printf("\nThere are less number of numbers in the file!\n");
            free(a);
            free(b);
            free(x);
            free(y);
            free(mass);
            free(answer);
            MPI_Finalize();
            return -6;
         }
    } else
        InputMatrix(n, k, a, b, first_row, max_rows);

    if (my_rank == 0) printf("Matrix A is:\n");
    OutputMatrix(n, n, m, a, x, my_rank, p);
    if (my_rank == 0) printf("\nVector b is:\n");
    OutputMatrix(1, n, m, b, x, my_rank, p);

    for (i = 0; i < n; ++i)
        mass[i] = i;

    MPI_Barrier (MPI_COMM_WORLD);
    t = MPI_Wtime();

    if (SolveSystem(n, mass, a, b, x, y, answer, my_rank, p) == -1){
        if (my_rank == 0) printf("\nThere are no solutions: matrix A is degenerate!\n");
        free(a);
        free(b);
        free(x);
        free(y);
        free(mass);
        free(answer);
        MPI_Finalize();
        return -7;
    }

    MPI_Barrier (MPI_COMM_WORLD);
    t = MPI_Wtime() - t;

    if(my_rank == 0){
        printf("\nThe answer is:\n");
        OutputAnswer(m, answer);
    }

    if (k == 0)
        fileInputMatrix(n, argv[4], a, b, x, my_rank, p);
    else
        InputMatrix(n, k, a, b, first_row, max_rows);

    error = SolutionError(n, a, b, x, y, answer, my_rank, p);
    if(my_rank == 0) printf("\n||Ax - b||/||b|| = %10.3e\n", error);

    if(my_rank == 0) printf("\nSolution accuracy is: %10.3e\n", SolutionAccuracy(n, answer));

    if (my_rank == 0) printf("\nTime spent = %lf\n\n", t);

    free(a);
    free(b);
    free(x);
    free(y);
    free(answer);
    free(mass);

    MPI_Finalize();

    return 0;
}

