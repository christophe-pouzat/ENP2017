#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s [--F_inf=real] [--Delta=real] [--beta=real] [--t_onset=real] ...\n" \
    "          ... [--t_0=real] [--delta_t=real] [--n_obs=int]\n\n"\
    "  --F_inf <positive real>: asymptotic value (default 800)\n"\
    "  --Delta <positive real>: signal jump amplitude (default 150)\n"\
    "  --beta <positive real>: inverse time constant (default 0.1)\n"\
    "  --t_onset <real>: the stimulus onset time (default 5)\n"	   \
    "  --t_0 <real>: observations start time (default 0)\n"\
    "  --delta_t <positive real>: sampling period (default 0.15)\n"\
    "  --n_obs <positive integer>: number of observations (default 168)\n\n";  
  double F_inf = 800;
  double Delta = 150;
  double beta = 0.1;
  double t_onset = 5;
  double t_0 = 0;
  double delta_t = 0.15;
  size_t n_obs = 168;
  {int opt;
    static struct option long_options[] = {
      {"F_inf",optional_argument,NULL,'f'},
      {"Delta",optional_argument,NULL,'d'},
      {"beta",optional_argument,NULL,'b'},
      {"t_onset",optional_argument,NULL,'o'},
      {"t_0",optional_argument,NULL,'s'},
      {"delta_t",optional_argument,NULL,'t'},
      {"n_obs",optional_argument,NULL,'n'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
  			    "hf:d:b:o:s:t:n",
  			    long_options,	
  			    &long_index)) != -1) {
      switch(opt) {
      case 'f':
      {
        F_inf = (double) atof(optarg);
        if (F_inf <= 0)
        {
  	fprintf(stderr,"The asymptotic value should be > 0.\n");
  	return -1;
        }
      }
      break;
      case 'd':
      {
        Delta = (double) atof(optarg);
        if (Delta <= 0)
        {
  	fprintf(stderr,"Delta should be > 0.\n");
  	return -1;
        }
      }
      break;
      case 'b':
      {
        beta = (double) atof(optarg);
        if (beta <= 0)
        {
  	fprintf(stderr,"beta should be > 0.\n");
  	return -1;
        } 
      }
      break;
      case 'o':
      {
        t_onset = (double) atof(optarg); 
      }
      break;
      case 's':
      {
        t_0 = (double) atof(optarg); 
      }
      break;
      case 't':
      {
        delta_t = (double) atof(optarg);
        if (delta_t <= 0)
        {
  	fprintf(stderr,"delta_t should be > 0.\n");
  	return -1;
        } 
      }
      break;
      case 'n':
      {
        n_obs = (size_t) atoi(optarg);
        if (n_obs <= 0)
        {
  	fprintf(stderr,"n_obs should be > 0.\n");
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
  gsl_rng * r;
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);
  for (size_t i=0; i<n_obs; i++) {
    double t = t_0 + i*delta_t;
    double ideal = F_inf;
    if (t >= t_onset) ideal += Delta*exp(-beta*(t-t_onset));
    double obs = gsl_ran_poisson(r, ideal);
    printf("%g %g %g\n",t,ideal,obs);
  }
  gsl_rng_free(r);
  return 0;
}
