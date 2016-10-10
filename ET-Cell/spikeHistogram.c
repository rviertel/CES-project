#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ETCell.h"

void histogram(FILE* fp, FILE* op);

#define SPIKE_POINT 0
#define PI 3.1415926535897932384626433832795028841971693993751

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

int main(int argc,char* argv[])
{
  FILE* fp = fopen("data/1hz50/ET_2_150_50.dat","r");
  FILE* op = fopen("1histogram_2_150_50.dat","w");

  histogram(fp,op);

  fclose(fp);
  fclose(op);
  return 0;
}

void histogram(FILE* fp, FILE* op)
{
  //input metrics
  double t, V, nK, hNaP, mH, mCaT, hCaT, wBK, Ca, nMystery, tau;

  int i;

  //flags
  int scan_status;

  //auxiliary quantities
  double prev_V;
  double spikeTimes[1000];
  int spike_count = 0;

  scan_status = fscanf (fp," %le %le %le %le %le %le %le %le %le %le %le\n", &t, &V, &nK, &hNaP, &mH, &mCaT, &hCaT, &wBK, &Ca, &nMystery, &tau);

  prev_V = V;

  while(scan_status != EOF)
  {
    scan_status = fscanf (fp," %le %le %le %le %le %le %le %le %le %le %le\n", &t, &V, &nK, &hNaP, &mH, &mCaT, &hCaT, &wBK, &Ca, &nMystery, &tau);

      //spike occurs
      if(prev_V < SPIKE_POINT && V > SPIKE_POINT && t > 1000)
      {
        spikeTimes[spike_count] = t;
        spike_count++;
      }

    prev_V = V;
  } // while(scan_status != EOF)

  printf("spike_count: %d\n", spike_count);

  for(i=0;i<spike_count;i++)
  {
    fprintf(op,"%.17e\n", spikeTimes[i]);
  }


} // histogram()
