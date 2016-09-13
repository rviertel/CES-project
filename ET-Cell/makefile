executables: ET sample periodic

periodic: periodic_main.c ET.c processTrace.c parameters.h
		cc -o periodic periodic_main.c ET.c processTrace.c -lgsl -lgslcblas -lm -fopenmp

ET: ET_main.c ET.c processTrace.c parameters.h
	cc -o ET ET_main.c ET.c processTrace.c -lgsl -lgslcblas -lm -fopenmp

sample:sample_main.c ET.c processTrace.c parameters.h
	cc -o sample sample_main.c ET.c processTrace.c -lgsl -lgslcblas -lm -fopenmp

.PHONY: clean
clean:
	rm ET sample periodic metrics.dat *.o
