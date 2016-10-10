objects = sample.o abscissa.o ML.o readTrace.o

uq: $(objects)
	mpicc -o uq $(objects) -lgsl -lgslcblas -lm -fopenmp

ML.o: ML.c
	cc -c ML.c -lgsl -lgslcblas -lm
readTrace.o: readTrace.c
	cc -c readTrace.c
sample.o: sample.c
	mpicc -c sample.c -fopenmp
abscissa.o: abscissa.c
	cc -c abscissa.c -lm

.PHONY: clean
clean:
	rm uq $(objects)
