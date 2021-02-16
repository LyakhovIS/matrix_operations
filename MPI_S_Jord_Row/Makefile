CFLAGS= -O3 -mfpmath=sse -fstack-protector-all -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Wcast-qual -Wno-suggest-attribute=format

all: program

program: main.o matrix.o task.o func.h
	mpicc $(CFLAGS) main.o matrix.o task.o func.h -o ./a.out

main.o: main.c
	mpicc -c $(CFLAGS)  main.c

matrix.o: matrix.c
	mpicc -c $(CFLAGS)  matrix.c

task.o: task.c
	mpicc -c $(CFLAGS)  task.c

clean:
	rm -rf *.o program
