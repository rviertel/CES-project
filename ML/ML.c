#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int func (double t, const double y[], double dy[],
      void *params) {

  double V = y[0];
  double w = y[1];
  double current = *(double*)params;

  double VK = -84.0, Vl = -60.0, VCa = 120.0;
  double gk = 8.0, gl = 2.0, c = 20.0;
  double v1 = -1.2, v2 = 18.0;
  double v3 = 12.0, v4 = 17.4, phi = 2.0/30.0, gca = 4.0;
  // double v3 = 2.0, v4 = 30.0, phi = .04, gca = 4.4;

  double minf = 0.5*(1.0+tanh((V-v1)/v2));
  double winf = 0.5*(1.0+tanh((V-v3)/v4));
  double tauw = 1.0/cosh((V-v3)/(2.0*v4));

  dy[0] = (current-gca*minf*(V-VCa)-gk*w*(V-VK)-gl*(V-Vl))/c;
  dy[1] = phi*(winf-w)/tauw;

  return GSL_SUCCESS;
}

int main(int argc, char* argv[]) {
  double input;
  sscanf(argv[1],"%lf",&input);

  gsl_odeiv2_system sys = {func, NULL, 2, &input};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				  1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 10000.0; //time in ms
  double y[2] = { -59.474, 0.00027 };
  FILE* fp = fopen("ML.dat","w");

  for (i = 0; i <= 40000; i++)
    {
      double ti = i * t1 / 40000.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS) {
      	  fprintf (stderr,"error, return value=%d\n", status);
      	  break;
      	}


      // output to file. 17 significant digits for full double precision
      fprintf (fp," %.17e %.17e %.17e\n", t, y[0], y[1]);
    }

  fclose(fp);
  gsl_odeiv2_driver_free (d);
  return 0;
}
