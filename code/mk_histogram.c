#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_histogram.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s --sample-size=integer [--min=real] [--max=real] ...\n"
    "          ... [--n_bins=integer] [--unormalize]\n\n"
    "  --sample-size <positive integer>: the sample size\n"									
    "  --min <double>: the left boundary of the histogram (default: min value of the data)\n"
    "  --max <double>: the right boundary of the histogram (default: max value of the data)\n"
    "  --n_bins <integer>: the number of bins (default: max-min/sample size/10)\n"
    "  --unormalize: if set the histogram is not normalized\n\n"
    " The program reads data from the stdin builts the histogram and prints it\n"
    " to the stdout.\n\n";
  size_t sample_size = 0;
  size_t n_bins;
  double min,max;
  int normalize=1,min_set=0,max_set=0,bins_set=0;
  {int opt;
    static struct option long_options[] = {
      {"sample_size",required_argument,NULL,'s'},
      {"min",optional_argument,NULL,'d'},
      {"max",optional_argument,NULL,'m'},
      {"n_bins",optional_argument,NULL,'b'},
      {"unormalize",no_argument,NULL,'u'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hub:s:d:m:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'u':
      {
        normalize = 0;
      }
      break;
      case 'b':
      {
        n_bins = (size_t) atoi(optarg);
        bins_set = 1;
      }
      break;
      case 's':
      {
        sample_size = (size_t) atoi(optarg);
      }
      break;
      case 'd':
      {
        min = (double) atof(optarg);
        min_set = 1;
      }
      break;
      case 'm':
      {
        max = (double) atof(optarg);
        max_set = 1;
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
    if (min_set & max_set & min>max) {
      fprintf(stderr,"min > max! Not allowed.\n");
      return -1;
    }
  }
  double data[sample_size];
  size_t data_idx=0;
  double x;
  while (fscanf(stdin, "%lg", &x) == 1 & \
         data_idx < sample_size) {
    data[data_idx] = x;
    data_idx++;
  }
  if (!min_set | !max_set) {
    gsl_vector_view data_v = gsl_vector_view_array (data, sample_size);
    double xmin,xmax;
    gsl_vector_minmax(&data_v.vector,&xmin,&xmax);
    if (!min_set)
      min = xmin;
    if (!max_set)
      max = xmax;
  }
  if (!bins_set)
    n_bins = (size_t) ceil((max-min)*sample_size/10);
  gsl_histogram *h = gsl_histogram_alloc (n_bins);
  gsl_histogram_set_ranges_uniform (h, min, max);
  for (size_t i=0; i<sample_size; i++)
    gsl_histogram_increment(h, data[i]);
  if (normalize)
    gsl_histogram_scale (h,(double)n_bins/(sample_size*(max-min)));
  gsl_histogram_fprintf (stdout, h, "%g", "%g");
  gsl_histogram_free (h);
  return 0;
}
