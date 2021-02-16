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

    for (i = 0; i < n; ++i)                                                         // сначала вычисляем норму матрицы: считаем норму первой строки
        maxnorm += fabs(a[i]);

    for (i = 1; i < max_rows; ++i){                                                 // сравниваем с нормами остальных строк, получаем норму подматрицы
        norm = 0;
        for (j = 0; j < n; ++j)
            norm += fabs(a[i * n + j]);
        if (norm > maxnorm)
            maxnorm = norm;
    }
    MPI_Allreduce (&maxnorm, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);         // считаем максимальную норму по подматрицам

    if (tmp < EPS){                                                                 // делим матрицу на нее, если необходимо
        for (i = 0; i < max_rows; ++i){
            for (j = 0; j < n; ++j)
                a[i * n + j] /= tmp;
            b[i] /= tmp;
        }
    }

    for (i = 0; i < n; ++i){
        if (i >= first_row && i <= last_row){                                       // сюда может войти только определенный процесс
            rank = my_rank;                                                         // поэтому проверка на my_rank не требуется
            maxi = fabs(a[(i - first_row) * n + i]);
            l = i;
            for (j = i + 1; j < n; ++j)                                             // поиск максимального элемента строки i
                if (fabs(a[(i - first_row) * n + j]) > maxi){
                    l = j;
                    maxi = fabs(a[(i - first_row) * n + j]);
                }

            for (k = 0; k < p; ++k){                                                // отсылаем всем процессам номер данного процесса
                if (k != my_rank)
                    MPI_Send(&my_rank, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
            }
            if (fabs(maxi) < EPS){                                                  // если встретилась нулевая строка, то завершаем данный процесс -1
                flag[0] = -1;                                                       // всем остальным процессам передаем в переменной flag[0] -1
                MPI_Bcast(&flag[0], 1, MPI_INT, rank, MPI_COMM_WORLD);
                return -1;
            }
            MPI_Bcast(&flag[0], 1, MPI_INT, rank, MPI_COMM_WORLD);

            for (k = 0; k < p; ++k)
                if (k != my_rank)
                    MPI_Send(&l, 1, MPI_INT, k, 0, MPI_COMM_WORLD);                 // посылаем всем остальным процессам номер столбца с макс элементом

            if (l != i){                                                            // если столбцы не совпадают, меняем столбец l на столбец i
                for (j = 0; j < max_rows; ++j){
                    tmp = a[j * n + i];
                    a[j * n + i] = a[j * n + l];
                    a[j * n + l] = tmp;
                }
                k = mass[i];                                                        // не забываем про массив индексов
                mass[i] = mass[l];
                mass[l] = k;
            }

            maxi = a[(i - first_row) * n + i];
            for (j = i, k = 0; j < n; ++j, ++k){                                    // делим всю i-ю строчку на [i,i] элемент
                a[(i - first_row) * n + j] /= maxi;
                x[k] = a[(i - first_row) * n + j];                                  // записываем все в буфер
            }
            b[i - first_row] /= maxi;                                               // не забываем про вектор
            x[k] = b[i - first_row];                                                // координату вектора записываем в последний элемент буффера
            for (k = 0; k < p; ++k)
                MPI_Send(&my_rank, 1, MPI_INT, k, 0, MPI_COMM_WORLD);               // посылаем всем процессам номер текущего процесса
            MPI_Bcast(x, n + 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);                  // посылаем всем процессам вектор, состоящий из всех элементов
                                                                                    // i строки, правее i-го элемента, последнее число - i элемент вектора

        } else{
            MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);// получаем номер искавшего макс элемент процесса
            MPI_Bcast(&flag[0], 1, MPI_INT, rank, MPI_COMM_WORLD);                  // получаем значение флага, если оно -1, то завершаем процесс -1
            if (flag[0] == -1)
                return -1;


            MPI_Recv(&l, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);  // получаем всеми процессами номер столбца с макс элементом

            if (l != i){                                                           // если не совпадают, меняем
                for (j = 0; j < max_rows; ++j){
                    tmp = a[j * n + i];
                    a[j * n + i] = a[j * n + l];
                    a[j * n + l] = tmp;
                }
                k = mass[i];                                                        // опять же индексы
                mass[i] = mass[l];
                mass[l] = k;
            }
        }

        MPI_Recv(&rank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);    // получаем номер процесса, из которого высылали буффер
        if (rank != my_rank)                                                        // получаем всеми остальными процессами буффер
            MPI_Bcast(x, n + 1, MPI_DOUBLE, rank, MPI_COMM_WORLD);

        for (j = first_row; j <= last_row; ++j){                                    // производим вычисления
            if (j != i){                                                            // x[n - i] - координата вектора, (последний элемент буффера)
                b[j - first_row] -= x[n - i] * a[(j - first_row) * n + i];          // опять же, j строка большой матрицы - j - first_row строка
                for (k = n - 1; k >= i; --k)                                        // матрицы процесса
                    a[(j - first_row) * n + k] -= a[(j - first_row) * n + i] * x[k - i];
            }
        }
    }

    if (my_rank == 0){                                                              // чтобы переставить в ответе индексы в соответствии с массивом
        for (i = 0; i < max_rows; ++i)                                              // нулевым процессом объединяем все части вектора в один большой
            y[i] = b[i];                                                            // сначала записываем часть нулевого процесса
        if (p > 1){
            for (k = 1; k < p; ++k){
                MPI_Recv(&first, 1, MPI_INT, k, 0, MPI_COMM_WORLD, &status);        // потом получаем части от других процессов вместе с
                MPI_Recv(&last, 1, MPI_INT, k, 1, MPI_COMM_WORLD, &status);         // first_row и last_row этих процессов
                MPI_Recv(x, n, MPI_DOUBLE, k, 2, MPI_COMM_WORLD, &status);

                for (i = first; i <= last; ++i)                                     // записываем в вектор эти части
                    y[i] = x[i - first];
            }
        }
    } else{
        MPI_Send(&first_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);                     // всеми процессами посылаем нулевому процессу части вектора
        MPI_Send(&last_row, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);                      // и first_row last_row
        MPI_Send(b, n, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }

    if (my_rank == 0)                                                               // в итоге имеем в нулевом процессе вектор у
        for (i = 0; i < n; ++i)                                                     // составленный из всех частей
            answer[mass[i]] = y[i];                                                 // переставляем в соответствии с индексами
                                                                                    // и получаем в 0 процессе вектор answer
    return 0;
}
