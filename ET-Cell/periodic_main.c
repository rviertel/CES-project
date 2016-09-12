#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


typedef struct pars
{
  // 0 - constant input
  // 1 - periodic input
  int type;

  // constant input
  double input;

  // periodic input
  double frequency;
  double amplitude;
  double duty_cycle;
} Pars;

void ET(Pars* pars, FILE* fp, FILE* cp);
// void ET(FILE* input, FILE* fp, FILE* cp);
void processTrace(FILE* fp, FILE* op, double input);



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

    Pars input_pars;
    input_pars.type = 1;
    input_pars.frequency = freq;
    input_pars.amplitude = amp;
    input_pars.duty_cycle = duty;

    ET(&input_pars,fp,cp);

    fclose(fp);
    FILE* op = fopen("metrics.dat","a+");
    fp = fopen("ET.dat","r");

    processTrace(fp, op, amp);

    fclose(fp);
    fclose(op);
  }
  else
    fprintf(stderr, "syntax: ./periodic frequency amplitude duty_cycle\n");


  return 0;
}
