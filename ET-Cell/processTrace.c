#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "ETCell.h"

#define SPIKE_POINT 0

void processTrace(FILE* fp, FILE* op, InputData* input_data)
{
  //input metrics
  double t, V, nK, hNaP, mH, mCaT, hCaT, wBK, Ca, nMystery, tau;

  //flags
  int scan_status;
  int burst_status;
  int bursting_has_begun = 0;
  int burst_start_count = 0;
  int burst_end_count = 0;

  //output metrics
  double burstDuration = 0;
  double burstFrequency = 0;
  double spikesPerBurst = 0;
  double mmp = DBL_MAX;
  double amp = 0;
  double average_peak_firing_rate = 0;
  double average_firing_rate = 0;
  double vector_strength = 0;
  double mean_phase = 0;

  //auxiliary quantities
  double prev_V, prev_nK;
  double burstStartTimes[1000];
  double burstEndTimes[1000];
  double spikeTimes[1000];
  double peak_firing_rate = 0;
  double phase = 0;
  double vector[] = {0,0};
  int spike_count = 0;
  int valid_spikes = 0;
  int data_point_count = 0;
  int buffer = 4;
  double BURST_POINT;

  //TODO: think about how to count when a burst starts. This could be a data integrity issue if not resolved before running the final simulation
  if((*input_data).type == 0)
    BURST_POINT = 0.25 + ((*input_data).input*0.0003);
  else
    BURST_POINT = 0.25;

  scan_status = fscanf (fp," %le %le %le %le %le %le %le %le %le %le %le\n", &t, &V, &nK, &hNaP, &mH, &mCaT, &hCaT, &wBK, &Ca, &nMystery, &tau);

  if(nK > BURST_POINT)
    burst_status = 2;
  else
    burst_status = 0;

  prev_V = V;
  prev_nK = nK;

  while(scan_status != EOF)
  {
    scan_status = fscanf (fp," %le %le %le %le %le %le %le %le %le %le %le\n", &t, &V, &nK, &hNaP, &mH, &mCaT, &hCaT, &wBK, &Ca, &nMystery, &tau);

    // burst begins
    if( (prev_nK <= BURST_POINT && nK > BURST_POINT && burst_status == 0) ||
        (prev_V < SPIKE_POINT && V > SPIKE_POINT && burst_status == 0) )
    {
      burst_status = 1;
      if(buffer <= 0)
      {
        bursting_has_begun = 1;
        burstStartTimes[burst_start_count] = t;
        burst_start_count++;

        if((*input_data).type == 1)//if periodic
        {
          double freq = (*input_data).frequency;
          double period = 1.0/freq;
          double t_sec = t/1000.0;
          phase = fmod(t_sec,period)/period;
        }
      }
    }
    // burst ends
    if(prev_nK > BURST_POINT && nK <= BURST_POINT && burst_status == 1)
    {
      burst_status = 0;
      buffer--;
      if(bursting_has_begun)
      {
        burstEndTimes[burst_end_count] = t;
        burst_end_count++;
        valid_spikes += spike_count;
        if((*input_data).type == 1)
        {
          mean_phase += phase;
          vector[0] += cos(phase);
          vector[1] += sin(phase);
        }
        int i;
        for(i = 0;i < spike_count - 1; i++)
        {
          double inst_rate = 1.0/(spikeTimes[i+1] - spikeTimes[i]);
          inst_rate = inst_rate*1000;
          if(inst_rate > peak_firing_rate)
            peak_firing_rate = inst_rate;
          average_firing_rate += inst_rate;
        }
        average_peak_firing_rate += peak_firing_rate;
        spike_count = 0;
      }
    }

    //bursting
    if(burst_status && bursting_has_begun)
    {
      //spike occurs
      if(prev_V < SPIKE_POINT && V > SPIKE_POINT)
      {
        spikeTimes[spike_count] = t;
        spike_count++;
      }
    }

    if(V < mmp)
      mmp = V;

    amp += V;
    data_point_count++;



    prev_V = V;
    prev_nK = nK;
  } // while(scan_status != EOF)


  int num_bursts = min(burst_start_count, burst_end_count);

  average_firing_rate = average_firing_rate/(valid_spikes - 1);
  average_peak_firing_rate = average_peak_firing_rate/num_bursts;
  mean_phase = mean_phase/num_bursts;
  vector[0] = vector[0]/num_bursts;
  vector[1] = vector[1]/num_bursts;
  vector_strength = sqrt(vector[0]*vector[0] + vector[1]*vector[1]);

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

  burstFrequency = (burstFrequency/(num_bursts - 1))*1000;

  spikesPerBurst = (double)valid_spikes/(double)num_bursts;

  amp = amp/data_point_count;

  if((*input_data).type == 0)
    fprintf(op,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", (*input_data).input, burstDuration,burstFrequency,spikesPerBurst, mmp, amp, average_peak_firing_rate, average_firing_rate);
  else if((*input_data).type == 1)
    fprintf(op,"%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n", (*input_data).frequency, (*input_data).amplitude, (*input_data).duty_cycle, burstDuration, burstFrequency, spikesPerBurst, mmp, amp, average_peak_firing_rate, average_firing_rate, vector_strength, mean_phase);


} // processTrace()
