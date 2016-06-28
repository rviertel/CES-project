#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

void ET(double input, FILE* fp, FILE* cp);
void processTrace(FILE* fp, FILE* op, double input);

int main(int argc,char* argv[]) {

  double input;
  FILE* fp = fopen("ET.dat","w");
  FILE* cp = fopen("current.dat","w");
  char* check;

  FILE *file;
  if ((file = fopen("metrics.dat", "r")))
  {
      fclose(file);
      remove("metrics.dat");
  }

  if(argc > 1)
  {
    input = strtod(argv[1], &check);
    ET(input,fp,cp);

    fclose(fp);
    FILE* op = fopen("metrics.dat","a+");
    fp = fopen("ET.dat","r");

    processTrace(fp, op, input);

    fclose(fp);
    fclose(op);
  }
  else
    fprintf(stderr, "syntax ./ET double\n");


  return 0;
}
