/*include guards */
#ifndef ML_H_INCLUDED
#define ML_H_INCLUDED

#include <stdio.h>

/* Prototypes for the functions */
void ML(void* pars, size_t num_pars, FILE* fp);

void readTrace(FILE* trace, FILE* spikes, FILE* metrics);

#endif
