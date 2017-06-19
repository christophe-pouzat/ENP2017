#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s [--n=int,int,int,...]\n\n"\
    "  --n <positive integers>: Poisson parameters separated by ','\n\n"
    " The program writes to the stdout the CDF of a scaled Poisson random variable\n"
    " that is, if X_n is Poisson with parameter n, the first column of the output\n"
    " contains integers from 0 to n and the second column contains the corresponding\n"
    " value of the CFD of Y_n = (X_n-n)/sqrt(n)\n\n"; 
  int *n_seq;
  size_t nb_n=0;
  {int opt;
    static struct option long_options[] = {
      {"n",required_argument,NULL,'n'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hn:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'n':
      {
        char *start = strdup(optarg); // duplicate optarg content
        char *running;
        running = start;
        char *token  = strsep(&running, ",");
        while (token != NULL) {
  	token = strsep (&running, ","); // split d_name at each ","
  	nb_n++;
        }
        free(start);
        n_seq = malloc(nb_n*sizeof(int));
        start = strdup(optarg); // duplicate optarg content again
        running = start;
        // Get the parameter of each Poisson distribution
        for (size_t i=0; i<nb_n; i++) {
  	token = strsep (&running, ",");
  	int para = atoi(token);
  	if (para <= 0) {
  	  fprintf(stderr,"A Poisson parameter should be > 0.\n");
  	  free(start);
  	  free(n_seq);
  	  return -1;
  	} else {
  	  n_seq[i] = para;
  	}
        }
        free(start);
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
  }
  printf("# %d Scaled Poisson distributions "
         "with parameters:\n#  %d",(int) nb_n, n_seq[0]);
  if (nb_n > 1) {
    for (size_t i=1; i<nb_n; i++)
      printf(", %d",n_seq[i]);
  }
  printf("\n#\n#\n");
  for (size_t i=0; i<nb_n; i++) {
    // Loop over the Poisson parameters
    int n = n_seq[i];
    double sn = sqrt(n);
    double inv_sn = 1.0/sn;
    unsigned int from = n-4*sn > 0 ? n-4*sn : 0;
    unsigned int to = n+5*sn;
    printf("# Scaled Poisson CDF with parameter: %d\n",n);
    for (unsigned int k=from; k <= to; k++) {
      printf("%g %g\n", ((double) k-n)*inv_sn, gsl_cdf_poisson_P(k,n)); 
    }
    printf("\n\n");  
  }
  free(n_seq);
  return 0;
}
