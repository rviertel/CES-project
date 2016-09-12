#include <stdio.h>
#include <omp.h>

void ET(double input, FILE* fp, FILE* cp);
void processTrace(FILE* fp, FILE* op, double input);

int main(int argc,char* argv[]) {

  int i;

  FILE *file;
  if ((file = fopen("metrics.dat", "r")))
  {
      fclose(file);
      remove("metrics.dat");
  }

  for(i = 21; i <= 30; i+=1)
  {
    char trace[100];
    char current[100];

    sprintf(trace,"ET.dat");
    sprintf(current,"current.dat");

    FILE* fp = fopen(trace, "w");
    FILE* cp = fopen(current, "w");

    ET(i*10.0,fp,cp);

    fclose(fp);

    fp = fopen(trace, "r");
    FILE* op = fopen("metrics.dat","a+");

    processTrace(fp, op, i*10.0 - 5);

    fclose(fp);
    fclose(op);
  }


  return 0;
}
