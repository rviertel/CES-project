#ifndef UQ_ETCELL2_H_INCLUDED
#define UQ_ETCELL2_H_INCLUDED

#include <stdio.h>


typedef struct inputdata
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
} InputData;

void uq_ET(InputData* pars, FILE* fp);
void uq_processTrace(FILE* fp, FILE* op);

#endif
