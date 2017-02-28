#include <math.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "ETCell.h"
#include "parameters.h"

#define TIME 5000
#define NUM_EQ 10

int vfield (double t, const double y[], double dy[], void *params);
double squarewave(double t, double freq, double amp, double duty_cycle);

/*
███████  ██████  ██    ██  █████  ██████  ███████ ██     ██  █████  ██    ██ ███████
██      ██    ██ ██    ██ ██   ██ ██   ██ ██      ██     ██ ██   ██ ██    ██ ██
███████ ██    ██ ██    ██ ███████ ██████  █████   ██  █  ██ ███████ ██    ██ █████
     ██ ██ ▄▄ ██ ██    ██ ██   ██ ██   ██ ██      ██ ███ ██ ██   ██  ██  ██  ██
███████  ██████   ██████  ██   ██ ██   ██ ███████  ███ ███  ██   ██   ████   ███████
            ▀▀
*/

double squarewave(double t, double freq, double amp, double duty_cycle)
{
  double t_sec = t/1000;
  double period = 1.0/freq;
  double phase = fmod(t_sec,period)/period;

  if(phase*100 < duty_cycle)
    return amp;
  else
    return 0;
}



/*
 ██████  ██████  ███████     ██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
██    ██ ██   ██ ██          ██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
██    ██ ██   ██ █████       ██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
██    ██ ██   ██ ██          ██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
 ██████  ██████  ███████     ██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/


int vfield (double t, const double y[], double dy[], void *params) {

  double V = y[0];
  double nK = y[1];
  double hNaP = y[2];
  double hH = y[3];
  double mLVA1 = y[4];
  double hLVA1 = y[5];
  double mLVA2 = y[6];
  double wBK = y[7];
  double Ca = y[8];
  double nNew = y[9];

  InputData* input_data = (InputData*)params;
  double Input, Iex = 0;


  if((*input_data).type == 0) //constant input
    Input = (*input_data).input;

  if((*input_data).type == 1) // periodic input
  {
    double freq, amp, duty_cycle;
    freq = (*input_data).frequency;
    amp = (*input_data).amplitude;
    duty_cycle = (*input_data).duty_cycle;
    Input = squarewave(t,freq,amp,duty_cycle);
    Iex = (*input_data).input;
  }


  // auxiliary quantities for BK channel
  double theta_wBK = -32.0 + 59.2*exp(-90.0*Ca) + 96.7*exp(-470.0*Ca);
  double p = 2.9 + 6.3*exp(-360*Ca);
  // double p = 45.0/(1.0+exp(-(Ca - .05)/.007)); // value used when modulating SPB with BK
  double s = -25.3 + 107.5*exp(-120.0*Ca);
  double f = 1.0/(10.0*(exp(-(V+100.0-s)/63.6)+exp((-150.0+(V+100.0-s))/63.6))) - 5.2;


  // infinity curves
  double mNa_inf = 1.0/(1.0+exp((V-theta_mNa)/sigma_mNa));
  double nK_inf = 1.0/(1.0+exp((V-theta_nK)/sigma_nK));
  double hNaP_inf = 1.0/(1.0+exp((V-theta_hNaP)/sigma_hNaP));
  double hH_inf = 1.0/(1.0+exp((V-theta_hH)/sigma_hH));
  double mLVA1_inf = 1.0/(1.0+exp((V-theta_mLVA1)/sigma_mLVA1));
  double hLVA1_inf = 1.0/(1.0+exp((V-theta_hLVA1)/sigma_hLVA1));

  double alpha = 1.11/(1.0+exp(-0.058*(V+10)));
  double beta = (0.02*(V+23.9))/(1.0-exp((V+23.9)/5));
  double mLVA2_inf = alpha/(alpha+beta);

  double mNaP_inf = 1.0/(1.0+exp((V-theta_mNaP)/sigma_mNaP));
  double mHVA_inf = 1.0/(1.0+exp((V-theta_mHVA)/sigma_mHVA));
  double wBK_inf = 1.0/(1.0+exp((V-theta_wBK)/sigma_wBK));
  double mNew_inf = 1.0/(1.0+exp((V-theta_mNew)/sigma_mNew));
  double nNew_inf = 1.0/(1.0+exp((V-theta_nNew)/sigma_nNew));

  //time constants
  double nNew_tau = 1000/(1.0+exp(-(V+35))) + 1000;
  double nK_tau = tau_nK/cosh((V-theta_nK)/(2.0*sigma_nK));
  double hNaP_tau = tau_hNaP/cosh((V-theta_hNaP)/(2.0*sigma_hNaP));
  double hH_tau = tau_hH_T*exp(delta_hH_T*(V-theta_hH_T)/sigma_hH_T)
                            / (1+exp((V-theta_hH_T)/sigma_hH_T));
  double mLVA1_tau = 1.0/(exp((V+68.03)/-27.68)+exp((V+39.08)/2.74)) + 2.29;
  double hLVA1_tau;

  if(V < -50)
  {
    hLVA1_tau = exp((V+770)/162.5);
  }
  else
  {
    hLVA1_tau = 37 + exp((V+27)/-6);
  }

  double mLVA2_tau = 1.0/(alpha+beta);
  double wBK_tau = -(p - 1.0)*(f - 0.2)/0.8 + wBK_base;


  // compute values for the currents
  double INa = gNa*(1.0-nK)*mNa_inf*mNa_inf*mNa_inf*(V-vNa);
  double INew = gNew*mNew_inf*nNew*(V - vK);
  double IK = gK*nK*nK*nK*nK*(V-vK);
  double IL = gL*(V-vL);
  double IH = gH*hH*(V-vH);
  double INaP = gNaP*mNaP_inf*hNaP*(V-vNa);
  double ILVA1 = gLVA1*mLVA1*mLVA1*mLVA1*hLVA1*(V-vCa);
  double ILVA2 = gLVA2*mLVA2*mLVA2*(V-vCa);
  double IHVA = gHVA*mHVA_inf*(V-vCa);
  double IBK = gBK*wBK*(V-vK);


  dy[0] = -(INa + IK + ILVA1 + ILVA2 + IH + INaP + IL + IHVA + IBK - Input - Iex)/C;
  dy[1] = (nK_inf-nK)/nK_tau;
  dy[2] = (hNaP_inf-hNaP)/hNaP_tau;
  dy[3] = (hH_inf-hH)/hH_tau;
  dy[4] = (mLVA1_inf-mLVA1)/mLVA1_tau;
  dy[5] = (hLVA1_inf-hLVA1)/hLVA1_tau;
  dy[6] = (mLVA2_inf-mLVA2)/mLVA2_tau;
  dy[7] = (wBK_inf - wBK)/wBK_tau;
  dy[8] = -Ca_buffer*10.0*(ILVA1 + ILVA2 + IHVA)/(Ca_z*F*d) + (Ca0 - Ca)/tau_Ca;
  dy[9] = (nNew_inf - nNew)/nNew_tau;

  return GSL_SUCCESS;
}



/*
███████  ██████  ██     ██    ██ ███████ ██████
██      ██    ██ ██     ██    ██ ██      ██   ██
███████ ██    ██ ██     ██    ██ █████   ██████
     ██ ██    ██ ██      ██  ██  ██      ██   ██
███████  ██████  ███████  ████   ███████ ██   ██
*/

void ET(InputData* input_data, FILE* fp, FILE* cp, double Iex) {

  if((*input_data).type == 1)
    (*input_data).input = Iex;

  double begin = omp_get_wtime();

  gsl_odeiv2_system sys = {vfield, NULL, NUM_EQ, input_data};

  gsl_odeiv2_driver * driver =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
          1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = TIME; //time in ms
   double y[NUM_EQ] = {-51.7045, 0.0531, 0.0604, 0.1720, 0.0460, 0.2084, 0.0460, 0.1292, 0.0005, 0};

  for (i = 0; i <= (4*TIME); i++)
    {
      double ti = i * t1 / (4*TIME);

      int status = gsl_odeiv2_driver_apply (driver, &t, ti, y);

      if (status != GSL_SUCCESS) {
          fprintf (stderr,"error, return value=%d\n", status);
          break;
        }

      //compute currents
      double V = y[0];
      double nK = y[1];
      double hNaP = y[2];
      double hH = y[3];
      double mLVA1 = y[4];
      double hLVA1 = y[5];
      double mLVA2 = y[6];
      double wBK = y[7];
      double Ca = y[8];
      double nNew = y[9];
      double tau = 1000/(1.0+exp(-(V+35))) + 200;

      double mNew = 1.0/(1.0+exp((V-theta_nK)/sigma_nK));
      double mNa_inf = 1.0/(1.0+exp((V-theta_mNa)/sigma_mNa));
      double mHVA_inf = 1.0/(1.0+exp(-(V + 10.0)/6.5));
      double mNaP_inf = 1.0/(1.0+exp((V-theta_mNaP)/sigma_mNaP));


      double INa = gNa*(1.0-nK)*mNa_inf*mNa_inf*mNa_inf*(V-vNa); //2
      double IK = gK*nK*nK*nK*nK*(V-vK); //3
      double IL = gL*(V-vL); //4
      double IH = gH*hH*(V-vH); //5
      double INaP = gNaP*mNaP_inf*hNaP*(V-vNa); //6
      double ILVA1 = gLVA1*mLVA1*mLVA1*mLVA1*hLVA1*(V-vCa);
      double ILVA2 = gLVA2*mLVA2*mLVA2*(V-vCa); //7
      double IHVA = gHVA*mHVA_inf*(V-vCa); //8
      double IBK = gBK*wBK*(V-vK); //9
      double INew = gNew*mNew*nNew*(V - vK); //10

      // output to file. 17 significant digits for full double precision
      fprintf (fp," %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                  t, V, nK, hNaP, hH, mLVA1, hLVA1, mLVA2, wBK, Ca, nNew, tau);

      fprintf (cp," %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                  t, INa, IK, IL, IH, INaP, ILVA1, ILVA2, IHVA, IBK, INew);


    }

  gsl_odeiv2_driver_free (driver);

  double end = omp_get_wtime();

  fprintf(stdout, "%f\n", end-begin);
}
