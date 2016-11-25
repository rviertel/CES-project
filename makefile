objects = sample.o abscissa.o uq_ET2.o uq_processTrace2.o

uq: $(objects)
	mpicc -o uq $(objects) -L/uufs/chpc.utah.edu/sys/installdir/gsl/2.1/lib -Wl,-rpath=/uufs/chpc.utah.edu/sys/installdir/gsl/2.1 -lgsl -lgslcblas -lm -fopenmp

uq_ET.o: uq_ET2.c
	cc -c uq_ET2.c -L/uufs/chpc.utah.edu/sys/installdir/gsl/2.1/lib -Wl,-rpath=/uufs/chpc.utah.edu/sys/installdir/gsl/2.1  -lgsl -lgslcblas -lm
uq_processTrace.o: uq_processTrace2.c
	cc -c uq_processTrace2.c
sample.o: sample.c
	mpicc -c sample.c -fopenmp
abscissa.o: abscissa.c
	cc -c abscissa.c -lm

.PHONY: clean
clean:
	rm uq $(objects)
	rm -rf data
