#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <pthread.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "func.h"


using namespace std;


typedef struct {
    double *a;
    double *x;
    double *answer;
    int *mass;
    int n;
    int thread_num;
    int total_threads;
    int *flag;
    bool deg;
    double time;
} ARGS;

static double thread_time = 0;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *Solution(void *p_arg) {
	ARGS *arg = (ARGS*)p_arg;
	double t1;

    synchronize(arg->total_threads);
	t1 = get_full_time();
	if (SolveSystem(arg->n, arg->mass, arg->a, arg->x, arg->answer, arg->thread_num, arg->total_threads, arg->flag) == -1){
		arg->deg = true;}
	synchronize(arg->total_threads);

	t1 = get_full_time() - t1;
	arg->time = t1;
	thread_time += t1;
	synchronize(arg->total_threads);

	return 0;
}


int main(int argc, char *argv[]){
    pthread_t *threads;
    ARGS *args;

    int n, m, k, total_threads, i, p = 0;
    int *mass = NULL;
    int *deg = NULL;
    int *flag = NULL;
    char *name = NULL;
    double *a = NULL;
    double *b = NULL;
    double *x = NULL;
    double *answer = NULL;

    cout << endl;

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
    if (!(convent3 >> k)){
        cout << "Error with argument k!" << endl;
        return -1;
    }
    stringstream convent4(argv[4]);
    if (!(convent4 >> total_threads)){
        cout << "Error with argument total threads!" << endl;
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

    if (k != 0 && k != 1 && k != 2 && k != 3 && k != 4 && k != 5){
        cout << "Wrong argument k!" << endl;
        return -1;
    }

    cout << "Your arguments are:" << endl;
    cout << "n - " << n << endl;
    cout << "m - " << m << endl;
    cout << "k - " << k << endl;
    cout << "total threads - " << total_threads << endl;
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
    b = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    flag = (int*)malloc(sizeof(int));
    answer = (double*)malloc(n * sizeof(double));
    mass = (int*)malloc(n * sizeof(int));
    deg = (int*)malloc(total_threads * sizeof(int));
    threads = (pthread_t*)malloc(total_threads * sizeof(pthread_t));
    args = (ARGS*)malloc(total_threads * sizeof(ARGS));

    if (a == NULL || b == NULL || x == NULL || mass == NULL || answer == NULL || threads == NULL || args == 0 || deg == 0){
        cout << "Not enough memory!" << endl;
        free(a);
        free(b);
        free(x);
        free(mass);
        free(answer);
        free(threads);
        free(args);
        free(deg);
        free(flag);
        return -1;
    }

    for (i = 0; i < n; i++)
        mass[i] = i;
    flag[0] = 0;

    if (k == 0){
        if (fileInputMatrix(n, name, a) == -1){
            cout << "There are less number of numbers in the file!" << endl;
            free(a);
            free(b);
            free(x);
            free(mass);
            free(answer);
            free(threads);
            free(args);
            free(deg);
            free(flag);
            return -1;
        } else if(fileInputMatrix(n, name, a) == -2){
            cout << "File did not open!" << endl;
            free(a);
            free(b);
            free(x);
            free(mass);
            free(answer);
            free(threads);
            free(args);
            free(deg);
            free(flag);
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


    for (i = 0; i < total_threads; ++i){
        args[i].a = a;
        args[i].x = x;
        args[i].answer = answer;
        args[i].mass = mass;
        args[i].n = n;
        args[i].thread_num = i;
        args[i].total_threads = total_threads;
        args[i].deg = false;
        args[i].time = 0;
        args[i].flag = flag;
    }

    for (i = 0; i < total_threads; ++i){
		if (pthread_create(threads + i, 0, Solution, args + i)){
			cout << "Cannot create thread " << i << endl;
			if (a) free(a);
			if (b) free(b);
			if (x) free(x);
			if (answer) free(answer);
			if (mass) free(mass);
			if (threads) free(threads);
			if (args) free(args);
			if (deg) free(deg);
			if (flag) free(flag);
			return -1;
		}
    }

	for (i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], 0)){
			cout << "Cannot wait thread " << i << endl;
			if (a) free(a);
			if (b) free(b);
			if (x) free(x);
			if (answer) free(answer);
			if (mass) free(mass);
			if (threads) free(threads);
			if (args) free(args);
			if (deg) free(deg);
			if (flag) free(flag);
			return -1;
		}

	for (i = 0; i < total_threads; i++){
		if(args[i].deg){
			p = -1;
			break;
		}
	}


	if (p == -1){
	    cout << "There are no solutions: matrix A is degenerate!" << endl;
        free(threads);
        free(args);
        free(a);
        free(b);
        free(x);
        free(mass);
        free(answer);
        free(deg);
        free(flag);
        return -1;
	}

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

    printf("Total full time: %lf", args[total_threads-1].time);
    printf("\t total threads time: %lf", thread_time);
    printf("\t per thread: %lf\n", thread_time / double(total_threads));
    cout << endl;

    free(threads);
    free(args);
    free(a);
    free(b);
    free(x);
    free(mass);
    free(answer);
    free(deg);
    free(flag);

    return 0;
}
