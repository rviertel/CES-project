uq is a program to sample a system of differential equations over Guassian quadrature points in a specified parameter space.

The user must specify the number of parameters and parameter ranges, as well as the number of blocks and order of the quadrature following the example in the script. The user must also explicitly code a way for the main method to interact with the system of differential equations by writing a function call to the system. There is no API.

sampling data, and final metric data will be stored in a directory named "data".

There are 3 optional arguments (though at least one is necessary) to run uq.

-s -- Sample over the parameter space by integrating at each Guassian quadrature point.
-p -- Process the raw traces of the ODE system, to get data that will be used in the integration
-i -- Integrate to compute statistical moments of the processed data. Mean and variance will be computed for each metric over each block of parameter space.

To build, use the makefile. Note that openMP and MPI are required. You may need to change the flags in the makefile to include any necessary flags to compile the system of ODEs along with the solver. For example, in the example script, -lgsl -lgslcblas, and -lm are used.

To run, type

mpirun -N int uq [flags]

where int specifies the number of nodes that will be used to sample and integrate, and [flags] are one of the three flags specified above.
