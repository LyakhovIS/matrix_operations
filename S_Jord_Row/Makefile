CFLAGS= -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format 

all: program

program: main.o matrix.o task.o func.h
	g++ $(CFLAGS) main.o matrix.o task.o func.h -o ./a.out

main.o: main.cpp
	g++ -c $(CFLAGS)  main.cpp

matrix.o: matrix.cpp
	g++ -c $(CFLAGS)  matrix.cpp

task.o: task.cpp
	g++ -c $(CFLAGS)  task.cpp

clean:
	rm -rf *.o program
