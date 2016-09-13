#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "ETCell.h"

int main(int argc,char* argv[]) {

  double input;
  FILE* fp = fopen("ET.dat","w");
  FILE* cp = fopen("current.dat","w");
  char* check;

  FILE* file;
  if ((file = fopen("metrics.dat", "r")))
  {
      fclose(file);
      remove("metrics.dat");
  }

  if(argc > 1)
  {
    input = strtod(argv[1], &check);

    InputData input_data;
    input_data.type = 0;
    input_data.input = input;

    ET(&input_data,fp,cp);

    fclose(fp);
    FILE* op = fopen("metrics.dat","a+");
    fp = fopen("ET.dat","r");

    processTrace(fp, op, &input_data);

    fclose(fp);
    fclose(op);
  }
  else
    fprintf(stderr, "syntax: ./ET input\n");


  return 0;
}
