#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; })

void processTrace(FILE* fp, FILE* op, double input)
{
  //flags
  int scan_status;
  int burst_status;
  int bursting_has_begun = 0;
  int burst_start_count = 0;
  int burst_end_count = 0;

  double t, V, nK, hNaP, mH, mCaT, hCaT, wBK, Ca, nMystery, tau;
  double prev_V, prev_nK;
  double burstStartTimes[100];
  double burstEndTimes[100];

  //metrics
  double burstDuration = 0;
  double burstFrequency = 0;
  double spikesPerBurst = 0;
  double mmp = DBL_MAX;
  double amp = 0;

  //auxiliary quantities
  double spike_count = 0;
  double valid_spikes = 0;
  double data_point_count = 0;
  int buffer = 4;
  // int burstcount;

  scan_status = fscanf (fp," %le %le %le %le %le %le %le %le %le %le %le\n",
          &t, &V, &nK, &hNaP, &mH, &mCaT, &hCaT, &wBK, &Ca, &nMystery, &tau);

  if(nK > 0.23)
    burst_status = 2;
  else
    burst_status = 0;

  prev_V = V;
  prev_nK = nK;

  while(scan_status != EOF)
  {
    scan_status = fscanf (fp," %le %le %le %le %le %le %le %le %le %le %le\n",
            &t, &V, &nK, &hNaP, &mH, &mCaT, &hCaT, &wBK, &Ca, &nMystery, &tau);

    // burst begins
    if(prev_nK <= 0.23 && nK > 0.23 && buffer <= 0)
    {
      bursting_has_begun = 1;
      burst_status = 1;
      burstStartTimes[burst_start_count] = t;
      burst_start_count++;

    }
    // burst ends
    if(prev_nK > 0.23 && nK <= 0.23)
    {
      burst_status = 0;
      buffer--;
      if(bursting_has_begun)
      {
        burstEndTimes[burst_end_count] = t;
        burst_end_count++;
        valid_spikes += spike_count;
        spike_count = 0;
      }
    }

    //bursting
    if(burst_status && bursting_has_begun)
    {
      //spike occurs
      if(prev_V < 0 && V > 0)
      {
        spike_count++;
      }
    }

    if(V < mmp)
      mmp = V;

    if(V < -40)
    {
      amp += V;
      data_point_count++;
    }



    prev_V = V;
    prev_nK = nK;
  } // while(scan_status != EOF)

  int num_bursts = min(burst_start_count, burst_end_count);

  int i;

  for(i = 0; i < num_bursts; i ++)
  {
    burstDuration += burstEndTimes[i] - burstStartTimes[i];
  }

  burstDuration = burstDuration/num_bursts;

  for(i = 0; i < num_bursts - 1; i ++)
  {
    burstFrequency += 1/(burstStartTimes[i+1] - burstStartTimes[i]);
  }

  burstFrequency = (burstFrequency/num_bursts)*1000;

  spikesPerBurst = valid_spikes/num_bursts;

  amp = amp/data_point_count;

  fprintf(op,"%.17e %.17e %.17e %.17e %.17e %.17e\n", input, burstDuration, burstFrequency,
                              spikesPerBurst, mmp, amp);




} // processTrace(char* filename)
