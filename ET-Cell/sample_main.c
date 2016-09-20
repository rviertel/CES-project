#include <stdio.h>
#include <omp.h>
#include "ETCell.h"

int main(int argc,char* argv[]) {

  InputData input_data;

  int i;

  FILE *file;
  if ((file = fopen("metrics.dat", "r")))
  {
      fclose(file);
      remove("metrics.dat");
  }

  for(i = 1; i <=10 ; i++)
  {
    char trace[100];
    char current[100];

    sprintf(trace,"ET.dat");
    sprintf(current,"current.dat");

    FILE* fp = fopen(trace, "w");
    FILE* cp = fopen(current, "w");

    input_data.type = 1;
    input_data.frequency = i;
    input_data.amplitude = 20;
    input_data.duty_cycle = 15;

    ET(&input_data,fp,cp);

    fclose(fp);

    fp = fopen(trace, "r");
    FILE* op = fopen("metrics.dat","a+");

    processTrace(fp, op, &input_data);

    fclose(fp);
    fclose(op);
  }


  return 0;
}
