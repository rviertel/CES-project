#include <stdio.h>
#include <omp.h>
#include "ETCell.h"

#define DUTY_CYCLE 70
#define MAX_i 25
#define MAX_j 10
#define MAX_k 5

int main(int argc,char* argv[]) {

  InputData input_data;

  int i,j,k;
  double ex_in[] = {-15,14,67,137,205};


  FILE *file;
  if ((file = fopen("data/metrics.dat", "r")))
  {
      fclose(file);
      remove("data/metrics.dat");
  }

  for(k = 1;k <= MAX_k; k++)
  {
    for(j = 1; j <= MAX_j; j++)
    {
      for(i = 0; i <= MAX_i ; i++)
      {
        char trace[100];
        char current[100];
        char metrics[100];

        int simulate = 0;

        sprintf(trace,"data/%dhz%d/ET_%d_%d.dat",k,DUTY_CYCLE,j,i*10);
        sprintf(current,"data/%dhz%d/current_%d_%d.dat",k,DUTY_CYCLE,j,i*10);
        sprintf(metrics,"data/%dhz%d/metrics.dat",k,DUTY_CYCLE);

        if (!(file = fopen(trace, "r")))
        {
          simulate = 1;
        }
        else
          fclose(file);


        input_data.type = 1;
        input_data.frequency = j;
        input_data.amplitude = 10*i;
        input_data.duty_cycle = DUTY_CYCLE;
        input_data.input = ex_in[k-1];

        if(simulate)
        {
          FILE* fp = fopen(trace, "w");
          FILE* cp = fopen(current, "w");
          ET(&input_data,fp,cp,ex_in[k-1]);
          fclose(fp);
          fclose(cp);
        }


        FILE* fp = fopen(trace, "r");
        FILE* op = fopen(metrics,"a+");

        processTrace(fp, op, &input_data);

        fclose(fp);
        fclose(op);

        printf("%d of %d iterations. frequency: %d amplitude: %d\n", (j-1)*(MAX_i+1) + i + 1 + (k-1)*(MAX_i+1)*MAX_j,(MAX_i+1)*MAX_j*MAX_k,j,i*10);
      }
    }
  }


  return 0;
}
