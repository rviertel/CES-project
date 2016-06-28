#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include "ML.h"

#define NUM_PARS 2
#define NUM_METRICS 2
#define EXTENSION_LEN 35
#define FILENAME_LEN 15 + EXTENSION_LEN

int main(int argc, char *argv[]) {

/***************initialize MPI environment**************************/
  int rank;
  int comm_size;

  MPI_Init(&argc,&argv);
  MPI_Comm world_comm = MPI_COMM_WORLD;
  MPI_Comm_size(world_comm,&comm_size);
  MPI_Comm_rank(world_comm,&rank);

/******************process arguments********************************/
  int i,j,sample = 0,process = 0,integrate = 0;

  // set flags
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i],"-s") == 0)
    sample = 1;
    if (strcmp(argv[i],"-p") == 0)
    process = 1;
    if (strcmp(argv[i],"-i") == 0)
    integrate = 1;
  }

  if (sample == 0 && process == 0 && integrate == 0) {
    fprintf(stderr,"Error: flags not set. Aborting...\n");
    return -1;
  }
/*******************define parameter ranges************************/
  typedef struct Parameters {
    char name[20];
    double low;
    double high;
    int blocks;
    double block_len;
    int order;
    int M;
  } Parameter;

Parameter pars[NUM_PARS];

// define parameters;
strcpy(pars[0].name, "current");
pars[0].low = 0;
pars[0].high = 100;
pars[0].blocks = 5;
pars[0].order = 8;
pars[0].block_len = (pars[0].high - pars[0].low)/pars[0].blocks; // don't change
pars[0].M = ((pars[0].order+1)>>1); // don't change

strcpy(pars[1].name, "poop");
pars[1].low = 50;
pars[1].high = 100;
pars[1].blocks = 2;
pars[1].order = 4;
pars[1].block_len = (pars[1].high - pars[1].low)/pars[1].blocks; // don't change
pars[1].M = ((pars[1].order+1)>>1); // don't change

/*****************define element grid***************************/
// define element structure
int num_elements = 1;
for(i=0;i<NUM_PARS;i++)
  num_elements = num_elements * pars[i].blocks;

  typedef struct Elements {
    double boundaries[NUM_PARS*2];
  } Element;

Element* elements = (Element*)malloc(num_elements*sizeof(Element));
int** index = (int**)malloc(num_elements*sizeof(int*));
for(i=0;i<num_elements;i++)
  index[i] = (int*)malloc(NUM_PARS*sizeof(int));

int count[NUM_PARS];

for(i=0;i<NUM_PARS;i++)
  count[i] = 0;

// sequentialize the index
// followed example from
// http://stackoverflow.com/questions/4683539/variable-amount-of-nested-for-loops
int k, pos;
for(k=0;k<num_elements;k++) {
  for(j=0;j<NUM_PARS;j++)
    index[k][j] = count[j];
  pos = NUM_PARS - 1;
  count[pos]++;
  while (pos >= 0 && count[pos] >= pars[pos].blocks) {
    count[pos] = 0;
    pos--;
    if (pos >= 0)
      count[pos]++;
  }
}

// define element boundaries
#pragma omp parallel for private(j)
for(k=0;k<num_elements;k++)
  for(j=0;j<NUM_PARS;j++) {
    elements[k].boundaries[j*2] = pars[j].block_len*index[k][j] + pars[j].low;
    elements[k].boundaries[j*2+1] = pars[j].block_len*(index[k][j]+1) + pars[j].low;
  }


FILE* efp = fopen("data/element_key.dat","w");

//print element boundaries
if (rank == 0)
  for(k=0;k<num_elements;k++) {
    fprintf(efp,"%d ", k);
    for(j=0;j<NUM_PARS;j++)
      fprintf(efp,"%d ", index[k][j]);
    for(j=0;j<NUM_PARS*2;j++)
      fprintf(efp,"%f ", elements[k].boundaries[j]);
    fprintf(efp,"\n");
  }

/**********************Assign work to each MPI task*****************/
// followed example on http://wiki.ccs.tulane.edu/index.php5/Parallel_Loop_MPI
int start, end;

start = (num_elements / comm_size) * rank;
 if (num_elements % comm_size > rank){
   start += rank;
   end = start + (num_elements / comm_size) + 1;
 }else{
   start += num_elements % comm_size;
   end = start + (num_elements / comm_size);
 }

 /******************set up quadrature grid index *****************************/
 //find number of grid points in each element
 int num_points = 1;
 for(j=0;j<NUM_PARS;j++)
  num_points = num_points*pars[j].order;

 //sequentialize the grid index
 int** grid_index = (int**)malloc(num_points*sizeof(int*));
 for(i=0;i<num_points;i++)
   grid_index[i] = (int*)malloc(NUM_PARS*sizeof(int));
 for(j=0;j<NUM_PARS;j++)
  count[j] = 0;

 for(i=0;i<num_points;i++) {
   for(j=0;j<NUM_PARS;j++)
    grid_index[i][j] = count[j];
   pos = NUM_PARS - 1;
   count[pos]++;
   while (pos >= 0 && count[pos] >= pars[pos].order) {
     count[pos] = 0;
     pos--;
     if (pos >= 0)
     count[pos]++;
   }
 }

 /*****************define abscissa and weights for each parameter*************/
 double* absc[NUM_PARS];
 double* weights[NUM_PARS];

 for(j=0;j<NUM_PARS;j++) {
   absc[j] = (double*)malloc(pars[j].M*sizeof(double));
   weights[j] = (double*)malloc(pars[j].M*sizeof(double));
   abscissa(pars[j].order,absc[j],weights[j]);
 }
/******************************************************************************/

  for (k = start; k < end; k++) {

/***************************find volume of this element************************/
double volume = 1;
for(j=0;j<NUM_PARS;j++)
  volume = volume*(elements[k].boundaries[j*2+1] - elements[k].boundaries[j*2]);

/************************* define quadrature grid for this element************/
  double quad_grid[num_points][NUM_PARS];
  double A[NUM_PARS],B[NUM_PARS];

  // Define A and B for each parameter
  for(j=0;j<NUM_PARS;j++) {
    A[j] = 0.5*(elements[k].boundaries[j*2+1] - elements[k].boundaries[j*2]);
    B[j] = 0.5*(elements[k].boundaries[j*2+1] + elements[k].boundaries[j*2]);
  }

  //assign values for each parameter at each grid point
  #pragma omp parallel for private(j)
  for(i=0;i<num_points;i++) {
    for(j=0;j<NUM_PARS;j++) {
      double X; //abscissa
      if (grid_index[i][j] < pars[j].M)
        X = -absc[j][pars[j].M - 1 - grid_index[i][j]];
      else {
        X = absc[j][grid_index[i][j] - pars[j].M + (pars[j].order&1)];
      }
      quad_grid[i][j] = A[j]*X+B[j];
    }
  }

  char estr[40];
  sprintf(estr,"data/grid_key_eid_%d.dat",k);
  FILE* gfp = fopen(estr,"w");

  // print grid points to file
  if (rank == 0)
    for(i=0;i<num_points;i++) {
      fprintf(gfp,"%d ", i);
      for(j=0;j<NUM_PARS;j++)
        fprintf(gfp,"%d ", grid_index[i][j]);
      for(j=0;j<NUM_PARS;j++)
        fprintf(gfp,"%f ", quad_grid[i][j]);
      fprintf(gfp,"\n");
    }
/******************************************************************************/


/*****************sample, process, and integrate within each element***********/
  double mean[NUM_METRICS], variance[NUM_METRICS];
  double sum[NUM_METRICS] = {0}, sumsum[NUM_METRICS] = {0};

  #pragma omp parallel private(j)
  {
  double sum_private[NUM_METRICS] = {0};
  double sumsum_private[NUM_METRICS] = {0};
  #pragma omp for
  for(i=0;i<num_points;i++) {

      double value[NUM_METRICS];
      char str[EXTENSION_LEN];
      char tracestr[FILENAME_LEN], processedstr[FILENAME_LEN];
      char metricstr[FILENAME_LEN];

      FILE *trace, *processed, *metric;

      // name data files
      sprintf(str, "_eid_%d_gid_%d.dat", k,i);

      sprintf(tracestr,"data/trace");
      strcat(tracestr, str);

      sprintf(processedstr,"data/processed");
      strcat(processedstr, str);

      sprintf(metricstr,"data/metric");
      strcat(metricstr, str);

      // sample
      if (sample == 1) {
        trace = fopen(tracestr, "w");
        ML(&quad_grid[i][0],NUM_PARS,trace);
        fclose(trace);
      }

      // process
      if (process == 1) {
        trace = fopen(tracestr, "r");
        processed = fopen(processedstr, "w");
        metric = fopen(metricstr, "w");
        read(trace, processed, metric);
        fclose(trace);
        fclose(processed);
        fclose(metric);
      }

      // integrate to find statistical moments
      if (integrate == 1) {
        metric = fopen(metricstr, "r");

        for(j=0;j<NUM_METRICS;j++)
          fscanf(metric,"%le",&value[j]);

        double weight = 1;
        for(j=0;j<NUM_PARS;j++) {
          if (grid_index[i][j] < pars[j].M)
            weight = weight * weights[j][pars[j].M-1 - grid_index[i][j]];
          else
            weight = weight * weights[j][grid_index[i][j] - (pars[j].M)+(pars[j].order&1)];
        }

        // sum over points that this thread is working on
        for(j=0;j<NUM_METRICS;j++) {
          sum_private[j] += value[j]*weight;
          sumsum_private[j] += value[j]*value[j]*weight;
        }

        fclose(metric);
      }

    }

    //add partial sums from each thread to get sum over all points
    #pragma omp critical
    {
        for(j=0;j<NUM_METRICS;j++) {
          sum[j] += sum_private[j];
          sumsum[j] += sumsum_private[j];
        }
    }
  }

  for(j=0;j<NUM_PARS;j++)
    for(i=0;i<NUM_METRICS;i++) {
      sum[i] = A[j]*sum[i];
      sumsum[i] = A[j]*sumsum[i];
    }

  for(i=0;i<NUM_METRICS;i++){
    mean[i] = sum[i]/volume;
    variance[i] = sumsum[i]/volume - mean[i]*mean[i];
  }

  for(i=0;i<NUM_METRICS;i++)
    printf("%.17le ", variance[i]);
  printf("\n");


}

  free(elements);

  for(i=0;i<num_elements;i++)
    free(index[i]);

  for(i=0;i<num_points;i++)
    free(grid_index[i]);

  for(j=0;j<NUM_PARS;j++) {
    free(absc[j]);
    free(weights[j]);
  }

  MPI_Finalize();
  return 0;
}
