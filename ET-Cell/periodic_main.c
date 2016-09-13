#include <stdlib.h>
#include <ctype.h>
#include "ETCell.h"


int main(int argc,char* argv[]) {

  double freq,amp,duty;
  FILE* fp = fopen("ET.dat","w");
  FILE* cp = fopen("current.dat","w");
  char* check;

  FILE* file;
  if ((file = fopen("metrics.dat", "r")))
  {
      fclose(file);
      remove("metrics.dat");
  }

  if(argc == 4)
  {
    freq = strtod(argv[1], &check);
    amp = strtod(argv[2], &check);
    duty = strtod(argv[3], &check);

    InputData input_data;
    input_data.type = 1;
    input_data.frequency = freq;
    input_data.amplitude = amp;
    input_data.duty_cycle = duty;

    ET(&input_data,fp,cp);

    fclose(fp);
    FILE* op = fopen("metrics.dat","a+");
    fp = fopen("ET.dat","r");

    processTrace(fp, op, &input_data);

    fclose(fp);
    fclose(op);
  }
  else
    fprintf(stderr, "syntax: ./periodic frequency amplitude duty_cycle\n");


  return 0;
}
