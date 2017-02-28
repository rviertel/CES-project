#include <stdio.h>
#include <omp.h>
#include "ETCell.h"

#define MAX_i 20
#define STEP 5

int main(int argc,char* argv[])
{

  InputData input_data;

  int i;

  FILE *file;
  if ((file = fopen("spb_data/metrics.dat", "r")))
  {
      fclose(file);
      remove("spb_data/metrics.dat");
  }

  for(i = 0; i <= MAX_i ; i++)
  {
    char trace[100];
    char current[100];
    char metrics[100];

    int simulate = 0;

    sprintf(trace,"spb_data/ET_%d.dat",i*STEP);
    sprintf(current,"spb_data/current_%d.dat",i*STEP);
    sprintf(metrics,"spb_data/metrics.dat");

    if (!(file = fopen(trace, "r")))
    {
      simulate = 1;
    }
    else
      fclose(file);


    input_data.type = 0;
    input_data.input = i*STEP;

    if(simulate)
    {
      FILE* fp = fopen(trace, "w");
      FILE* cp = fopen(current, "w");
      ET(&input_data,fp,cp,0);
      fclose(fp);
      fclose(cp);
    }


    FILE* fp = fopen(trace, "r");
    FILE* op = fopen(metrics,"a+");

    processTrace(fp, op, &input_data);

    fclose(fp);
    fclose(op);

    printf("%d of %d iterations. amplitude: %d\n", i+1, MAX_i, i*STEP);
  }

  return 0;
}
