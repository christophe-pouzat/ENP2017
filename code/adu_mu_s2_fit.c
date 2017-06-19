#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s --sample_size=integer\n\n"\
    "  --sample_size <positive integer>: the number of pixels times the number of exposures\n\n"
    " The program reads data from the stdin and performs a weighter linear fit.\n"
    " It prints to the stdout the results\n\n";
  size_t sample_size;
  {int opt;
    static struct option long_options[] = {
      {"sample_size",required_argument,NULL,'s'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hs:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 's':
      {
        sample_size = (size_t) atoi(optarg);
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
  }
  double mean[sample_size];
  double variance[sample_size];
  double w[sample_size];
  size_t data_idx=0;
  double mu,s2;
  while (fscanf(stdin, "%lg %lg", &mu, &s2) == 2 &	\
         data_idx < sample_size) {
    mean[data_idx] = mu;
    variance[data_idx] = s2;
    w[data_idx] = 0.5*99.0/s2/s2;
    data_idx++;
  }
  double c0, c1, cov00, cov01, cov11, chisq;
  
  gsl_fit_wlinear(mean, 1, w, 1, variance, 1, sample_size, 
  		&c0, &c1, &cov00, &cov01, &cov11, 
  		&chisq);
  printf ("Best fit: Variance(ADU) = %g + %g Mean(ADU)\n", c0, c1);
  printf ("Covariance matrix:\n");
  printf ("[ %g, %g\n  %g, %g]\n", 
          cov00, cov01, cov01, cov11);
  printf ("Residual sum of squares = %g\n", chisq);
  double gain_se = sqrt(cov11)*gsl_cdf_tdist_Pinv(0.975,(double)sample_size-3.0);
  printf ("The estimated gain is: %g\n",c1);
  printf ("A 0.95 confidence interval for the gain is:\n (%g,%g)\n",
          c1-gain_se,c1+gain_se);
  s2 = c0/gsl_pow_2(c1);
  printf ("The estimated read-out variance is: %g\n",s2);
  double s2_se = sqrt(cov00/gsl_pow_2(c1)+cov11*gsl_pow_2(2*c0/gsl_pow_3(c1)));
  s2_se *= gsl_cdf_tdist_Pinv(0.975,(double)sample_size-3.0);
  printf ("A 0.95 confidence interval for the read-out variance is:\n (%g,%g)\n",
          s2-s2_se,s2+s2_se);
  return 0;
}
