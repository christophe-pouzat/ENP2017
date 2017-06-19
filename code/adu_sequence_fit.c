#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s [--nb_frames=integer] [--conf_level=real]\n\n"\
    "  --nb_frames <positive integer>: the number of frames (default at 100)\n"
    "  --conf_level <positive real>: confidence level (between 0 and 1, default 0.95)\n\n"
    " The program reads data from the stdin and performs a linear fit.\n"
    " It prints to the stdout data compatible with gnuplot with the\n"
    " fitting results as comments followed by 4 columns:\n"
    "   - the frame index\n"
    "   - the fitted ADU\n"
    "   - the lower bound of the fitted ADU\n"
    "   - the upper bound of the fitted ADU.\n\n";
  size_t n = 100;
  double cl = 0.95;
  {int opt;
    static struct option long_options[] = {
      {"nb_frames",optional_argument,NULL,'n'},
      {"conf_level",optional_argument,NULL,'c'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hn:c:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'n':
      {
        n = (size_t) atoi(optarg);
      }
      break;
      case 'c':
      {
        cl = (double) atof(optarg);
        if (cl <= 0 | cl >= 1)
        {
          fprintf(stderr,"A confidence interval should be between 0 and 1.\n");
          return -1;
        }
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
  }
  double adu[n];
  {
    char input[256];
    size_t i=0;
    while (fgets(input,sizeof(input),stdin) != NULL & i < n) {
      adu[i] = (double) atoi(input);
      i++;
    }
    if (i != n)
      fprintf(stderr,"Fewer ADU than planed.\n");
  }
  double x[n];
  for (size_t i=0; i<n; i++)
    x[i] = (double) i;
  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear (x, 1, adu, 1, n, &c0, &c1,
  		&cov00, &cov01, &cov11, 
  		&sumsq);
  printf ("# best fit: ADU = %g + %g F\n", c0, c1);
  printf ("# covariance matrix:\n");
  printf ("# [ %g, %g\n#   %g, %g]\n", 
  	cov00, cov01, cov01, cov11);
  printf ("# Residual sum of squares = %g\n", sumsq);
  double slope_se = sqrt(cov11)*gsl_cdf_tdist_Pinv(1.0-(1.0-cl)*0.5,(double)n-3.0);
  printf ("# A %g confidence interval for the slope is (%g,%g)\n",
  	cl,c1-slope_se,c1+slope_se);
  printf ("#\n");
  printf ("# Frame index  Fitted value  Lower bound  Upper bound\n"); 
  double se_f = gsl_cdf_gaussian_Pinv(1.0-(1-cl)*0.5,1.0);
  for (size_t i=0; i<n; i++) {
    double aduf, aduf_err;
    gsl_fit_linear_est(x[i],c0,c1,cov00, cov01,
  		     cov11,&aduf,&aduf_err);
    printf("%g %g %g %g\n",
  	 x[i],aduf,aduf-se_f*aduf_err,
  	 aduf+se_f*aduf_err);
  }  
  return 0;
}
