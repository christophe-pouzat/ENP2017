#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_statistics_double.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s [--file=string] [--gain=real] [--variance=real]\n\n"\
    "  --file <character string>: data file name (default CCD_calibration.hdf5)\n"
    "  --gain <positive real>: CCD gain (default 1)\n"
    "  --variance <positive or null real>: CCD gain (default 0)\n\n"
    " The program opens 'file', and for each dataset, for each pixel\n"
    " computes the mean ADU and the ADU variance. If 'gain' and 'variance'\n"
    " are diffferent from 1 and 0 repectively then variance stabilization is\n"
    " performed, that is, the program works with 2*sqrt(ADU/gain+variance)\n"
    " instead of working with the ADU directly.\n"
    " The program prints its results to the stdout with the mean on the first\n"
    " column and the variance on the second.\n\n";
  size_t frame = 0;
  char default_file_name[] = "CCD_calibration.hdf5";
  char *filename;
  double gain = 1.0;
  double read_out_var = 0.0;
  {int opt;
    static struct option long_options[] = {
      {"file",optional_argument,NULL,'f'},
      {"gain",optional_argument,NULL,'g'},
      {"variance",optional_argument,NULL,'v'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    int file_name_set = 0;
    while ((opt = getopt_long(argc,argv,
                              "hf:g:v:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'f':
      {
        filename = optarg;
        file_name_set = 1;
      }
      break;
      case 'g':
      {
        gain = (double) atof(optarg);
        if (gain <= 0) {
          fprintf(stderr,"gain must be > 0\n");
          return -1;
        }
      }
      break;
      case 'v':
      {
        read_out_var = (double) atof(optarg);
        if (read_out_var < 0) {
          fprintf(stderr,"variance must be >= 0\n");
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
    if (file_name_set == 0) filename =  default_file_name;
  }
  // Open FILE
  hid_t file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  // Get dimensions of 3D object contained in DSET
  size_t n_groups = 10;
  char *groups[] = {"10ms","20ms","30ms","40ms","50ms",
  		  "60ms","70ms","80ms","90ms","100ms"};
  char *groupname = groups[0];
  char dset[25];
  dset[0] = '/';
  dset[1] = '\0';
  strcat(dset,groupname);
  strcat(dset,"/stack");
  hsize_t dims[3];
  H5LTget_dataset_info(file_id,dset,dims,NULL,NULL);
  size_t nrow = (size_t) dims[0];
  size_t ncol = (size_t) dims[1];
  size_t nsamp = (size_t) dims[2];
  // Read dset in data
  int *data = malloc(nrow*ncol*nsamp*sizeof(int));
  for (size_t g_idx=0; g_idx<n_groups; g_idx++) {
    char *groupname = groups[g_idx];
    dset[0] = '/';
    dset[1] = '\0';
    strcat(dset,groupname);
    strcat(dset,"/stack");
    H5LTread_dataset_int(file_id,dset,data);
    for (size_t i=0; i<nrow; i++) {
      for (size_t j=0; j<ncol; j++) {
        double y[nsamp];
        if (gain == 1.0 & read_out_var == 0.0) {
  	for (size_t k=0; k<nsamp; k++)
  	  // No variance stabilization
  	  y[k] = (double) data[(i*ncol+j)*nsamp+k];
        } else {
  	for (size_t k=0; k<nsamp; k++)
  	  // Do variance stabilization
  	  y[k] = 2*sqrt(((double) data[(i*ncol+j)*nsamp+k])/gain+read_out_var);
        }
        double mean = gsl_stats_mean(y, 1, nsamp);
        double variance = gsl_stats_variance_m(y, 1, nsamp, mean);
        printf("%g %g\n",mean,variance);
      }
    }
  }
  free(data);
  H5Fclose (file_id);  
  return 0;
}
