executables: ET sample

ET: ET_main.c ET.c processTrace.c parameters.h
	cc -o ET ET_main.c ET.c processTrace.c -lgsl -lgslcblas -lm -fopenmp

sample:sample.c ET.c processTrace.c parameters.h
	cc -o sample sample.c ET.c processTrace.c -lgsl -lgslcblas -lm -fopenmp

.PHONY: clean
clean:
	rm ET metrics.dat
