#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "parameters.h"

#define TIME 5000.0
#define NUM_EQ 9

int vfield (double t, const double y[], double dy[], void *params);


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
  double mH = y[3];
  double mCaT = y[4];
  double hCaT = y[5];
  double wBK = y[6];
  double Ca = y[7];
  double nMystery = y[8];

  double* pars = (double*)params;

  // get input
  double Input = pars[0];


  // auxiliary quantities for BK channel
  double theta_wBK = -32.0 + 59.2*exp(-90.0*Ca) + 96.7*exp(-470.0*Ca);
  double p = 2.9 + 6.3*exp(-360*Ca);
  double s = -25.3 + 107.5*exp(-120.0*Ca);
  double f = 1.0/(10.0*(exp(-(V+100.0-s)/63.6)+exp((-150.0+(V+100.0-s))/63.6))) - 5.2;


  // infinity curves
  double mNa_inf = 1.0/(1.0+exp((V-theta_mNa)/sigma_mNa));
  double nK_inf = 1.0/(1.0+exp((V-theta_nK)/sigma_nK));
  double hNaP_inf = 1.0/(1.0+exp((V-theta_hNaP)/sigma_hNaP));
  double mH_inf = 1.0/(1.0+exp((V-theta_mH)/sigma_mH));
  double mCaT_inf = 1.0/(1.0+exp((V-theta_mCaT)/sigma_mCaT));
  double hCaT_inf = 1.0/(1.0+exp((V-theta_hCaT)/sigma_hCaT));
  double mNaP_inf = 1.0/(1.0+exp((V-theta_mNaP)/sigma_mNaP));
  double mHVA_inf = 1.0/(1.0+exp(-(V + 10.0)/6.5));
  double wBK_inf = 1.0/(1.0+exp((theta_wBK - V)/15.6));

  // mystery equations
  double mMystery = 1.0/(1.0+exp((V+40)/-2));
  double nMystery_inf = 1.0/(1.0+exp((V+30)/-2));
  double nMystery_tau = 1000/(1.0+exp(-(V+35))) + 1000;

  //time constants
  double nK_tau = tau_nK/cosh((V-theta_nK)/(2.0*sigma_nK));
  double hNaP_tau = tau_hNaP/cosh((V-theta_hNaP)/(2.0*sigma_hNaP));
  double mH_tau = tau_mH_T*exp(delta_mH_T*(V-theta_mH_T)/sigma_mH_T)
                            / (1+exp((V-theta_mH_T)/sigma_mH_T));
  double mCaT_tau = tau_mCaT/cosh((V-theta_mCaT)/(2.0*sigma_mCaT));
  double hCaT_tau = tau_hCaT/cosh((V-theta_hCaT)/(2.0*sigma_hCaT));
  double wBK_tau = -(p - 1.0)*(f - 0.2)/0.8 + wBK_base;


  // compute values for the currents
  double INa = gNa*(1.0-nK)*mNa_inf*mNa_inf*mNa_inf*(V-vNa);
  double Imystery = gMystery*mMystery*nMystery*(V - vK);
  double IK = gK*nK*nK*nK*nK*(V-vK);
  double IL = gL*(V-vL);
  double IH = gH*mH*(V-vH);
  double INaP = gNaP*mNaP_inf*hNaP*(V-vNa);
  double ICaT = gCaT*mCaT*mCaT*hCaT*(V-vCa);
  double IHVA = gHVA*mHVA_inf*(V-vCa);
  double IBK = gBK*wBK*(V-vK);


  dy[0] = -(INa + IK + ICaT + IH + INaP + IL + IHVA + IBK + Imystery - Input)/C;
  dy[1] = (nK_inf-nK)/nK_tau;
  dy[2] = (hNaP_inf-hNaP)/hNaP_tau;
  dy[3] = (mH_inf-mH)/mH_tau;
  dy[4] = (mCaT_inf-mCaT)/mCaT_tau;
  dy[5] = (hCaT_inf-hCaT)/hCaT_tau;
  dy[6] = (wBK_inf - wBK)/wBK_tau;
  dy[7] = -Ca_buffer*10.0*(ICaT + IHVA)/(Ca_z*F*d) + (Ca0 - Ca)/tau_Ca;
  dy[8] = (nMystery_inf - nMystery)/nMystery_tau;

  return GSL_SUCCESS;
}

// int
// jac (double t, const double y[], double *dfdy,
//      double dfdt[], void *params)
// {
//   (void)(t); /* avoid unused parameter warning */
//   double mu = *(double *)params;
//   gsl_matrix_view dfdy_mat
//     = gsl_matrix_view_array (dfdy, 8, 8);
//   gsl_matrix * m = &dfdy_mat.matrix;
//   gsl_matrix_set (m, 0, 0, 0.0);
//   gsl_matrix_set (m, 0, 1, 1.0);
//   gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
//   gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
//   dfdt[0] = 0.0;
//   dfdt[1] = 0.0;
//   return GSL_SUCCESS;
// }

/*
███████  ██████  ██     ██    ██ ███████ ██████
██      ██    ██ ██     ██    ██ ██      ██   ██
███████ ██    ██ ██     ██    ██ █████   ██████
     ██ ██    ██ ██      ██  ██  ██      ██   ██
███████  ██████  ███████  ████   ███████ ██   ██
*/


void ET(double input, FILE* fp, FILE* cp) {

  double begin = omp_get_wtime();

  gsl_odeiv2_system sys = {vfield, NULL, NUM_EQ, &input};

  gsl_odeiv2_driver * driver =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk2,
          1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = TIME; //time in ms
   double y[NUM_EQ] = {-51.7045, 0.0531, 0.0604, 0.1720, 0.0460, 0.2084, 0.1292, 0.0005, 0};

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
      double mH = y[3];
      double mCaT = y[4];
      double hCaT = y[5];
      double wBK = y[6];
      double Ca = y[7];
      double nMystery = y[8];
      double tau = 1000/(1.0+exp(-(V+35))) + 200;

      double mMystery = 1.0/(1.0+exp((V-theta_nK)/sigma_nK));
      double mNa_inf = 1.0/(1.0+exp((V-theta_mNa)/sigma_mNa));
      double mHVA_inf = 1.0/(1.0+exp(-(V + 10.0)/6.5));
      double mNaP_inf = 1.0/(1.0+exp((V-theta_mNaP)/sigma_mNaP));


      double INa = gNa*(1.0-nK)*mNa_inf*mNa_inf*mNa_inf*(V-vNa); //2
      double IK = gK*nK*nK*nK*nK*(V-vK); //3
      double IL = gL*(V-vL); //4
      double IH = gH*mH*(V-vH); //5
      double INaP = gNaP*mNaP_inf*hNaP*(V-vNa); //6
      double ICaT = gCaT*mCaT*mCaT*hCaT*(V-vCa); //7
      double IHVA = gHVA*mHVA_inf*(V-vCa); //8
      double IBK = gBK*wBK*(V-vK); //9
      double Imystery = gMystery*mMystery*nMystery*(V - vK); //10

      // output to file. 17 significant digits for full double precision
      fprintf (fp," %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                  t, V, nK, hNaP, mH, mCaT, hCaT, wBK, Ca, nMystery, tau);

      fprintf (cp," %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                  t, INa, IK, IL, IH, INaP, ICaT, IHVA, IBK, Imystery);


    }

  gsl_odeiv2_driver_free (driver);

  // const gsl_odeiv_step_type * T
  //   = gsl_odeiv_step_rk8pd;
  //
  // gsl_odeiv_step * s
  //   = gsl_odeiv_step_alloc (T, NUM_EQ);
  // gsl_odeiv_control * c
  //   = gsl_odeiv_control_y_new (1e-6, 0.0);
  // gsl_odeiv_evolve * e
  //   = gsl_odeiv_evolve_alloc (NUM_EQ);
  //
  // double input;
  //
  // sscanf(argv[1],"%lf",&input);
  //
  // gsl_odeiv_system sys = {vfield, NULL, NUM_EQ, &input};
  //
  // int i;
  // int status;
  // double t = 0.0, t1 = 1000.0; //time in ms
  // double h = 1e-6;
  // double y[NUM_EQ] = {-51.7045, 0.0531, 0.0604, 0.1720, 0.0460, 0.2084, 0.1292, 0.0005};
  // FILE* fp = fopen("ET.dat","w");
  //
  // for (i = 0; i <= 12000; i++)
  //   {
  //     double ti = i * t1 / 12000.0;
  //     while (t < ti)
  //       {
  //         status = gsl_odeiv_evolve_apply (e, c, s,
  //                                              &sys,
  //                                              &t, t1,
  //                                              &h, y);
  //
  //         if (status != GSL_SUCCESS)
  //             break;
  //
  //
  //       }
  //       // output to file. 17 significant digits for full double precision
  //       fprintf(fp, "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n"
  //       ,t, y[0], y[1],y[2], y[3], y[4], y[5], y[6], y[7]);
  //   }
  //
  //
  // gsl_odeiv_evolve_free (e);
  // gsl_odeiv_control_free (c);
  // gsl_odeiv_step_free (s);

  double end = omp_get_wtime();

  fprintf(stdout, "%f\n", end-begin);
}
