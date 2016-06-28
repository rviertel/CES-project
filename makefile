objects = ML.o read.o sample.o abscissa.o

test: $(objects)
	mpicc -o test $(objects) -lgsl -lgslcblas -lm -fopenmp

ML.o: ML.c
	cc -c ML.c -lgsl -lgslcblas -lm
read.o:
sample.o: sample.c
	mpicc -c sample.c -fopenmp
abscissa.o: abscissa.c
	cc -c abscissa.c -lm

.PHONY: clean
clean:
	rm test $(objects)
