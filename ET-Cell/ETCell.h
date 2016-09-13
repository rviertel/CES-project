#ifndef ETCELL_H_INCLUDED
#define ETCELL_H_INCLUDED

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

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

void ET(InputData* pars, FILE* fp, FILE* cp);
void processTrace(FILE* fp, FILE* op, InputData* input_data);

#endif
