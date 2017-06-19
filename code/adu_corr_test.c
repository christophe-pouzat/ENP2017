#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <gsl/gsl_statistics_double.h>
#include <hdf5.h>
#include <hdf5_hl.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s [--group=string] [--file=string]\n\n"\
    "  --group <character string>: the name of the group (eg, 10ms the default)\n"
    "  --file <character string>: data file name (default CCD_calibration.hdf5)\n\n"
    " The program opens 'file', gets to dataset 'group/stack', gets the correlation\n"
    " coefficient of every pair of nearest neighbor ADU sequences and prints them to\n"
    " the stdout.\n"
    " Under the null hypothesis these correlation coefficients follow a Gaussian\n"
    " distribution centered on 0 with a variance of 1/100.\n\n";
  size_t frame = 0;
  char default_file_name[] = "CCD_calibration.hdf5";
  char *filename;
  char dset[25];
  dset[0] = '/';
  dset[1] = '\0';
  {int opt;
    size_t n_groups = 10;
    char *groups[] = {"10ms","20ms","30ms","40ms","50ms",
                      "60ms","70ms","80ms","90ms","100ms"};
    char *groupname;
    static struct option long_options[] = {
      {"group",optional_argument,NULL,'g'},
      {"file",optional_argument,NULL,'f'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    int file_name_set = 0, group_name_set = 0;
    while ((opt = getopt_long(argc,argv,
                              "hf:g:",
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
        groupname = optarg;
        size_t j;
        for (j=0; j<n_groups; j++)
          if (strcmp(groupname,groups[j])==0)
            break;
        if (j == n_groups) {
          fprintf(stderr,"Unknown group\n");
          return -1;
        }
        group_name_set = 1;
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
    if (group_name_set == 0) groupname =  groups[0];
    if (file_name_set == 0) filename =  default_file_name;
    strcat(dset,groupname);
    strcat(dset,"/stack");
  }
  // Open FILE
  hid_t file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  // Get dimensions of 3D object contained in DSET
  hsize_t dims[3];
  H5LTget_dataset_info(file_id,dset,dims,NULL,NULL);
  size_t nrow = (size_t) dims[0];
  size_t ncol = (size_t) dims[1];
  size_t nsamp = (size_t) dims[2];
  // Read dset in data
  int *data = malloc(nrow*ncol*nsamp*sizeof(int));
  H5LTread_dataset_int(file_id,dset,data);
  double corr;
  for (size_t i=1; i<(nrow-1); i+=2) {
    for (size_t j=1; j<(ncol-1); j+=2) {
      double adu_ref[nsamp];
      for (size_t idx=0; idx<nsamp; idx++)
        adu_ref[idx] = (double) data[(i*ncol+j)*nsamp+idx];
      double adu_test[nsamp];
      for (size_t idx=0; idx<nsamp; idx++) 
        adu_test[idx] = (double) data[((i-1)*ncol+j)*nsamp+idx];
      // look at pixel one row above
      corr = gsl_stats_correlation(adu_ref, 1, adu_test, 1, nsamp);
      printf("%g\n",corr);
      for (size_t idx=0; idx<nsamp; idx++)
        adu_test[idx] = (double) data[((i+1)*ncol+j)*nsamp+idx];
      // look at pixel one row below
      corr = gsl_stats_correlation(adu_ref, 1, adu_test, 1, nsamp);
      printf("%g\n",corr);
      for (size_t idx=0; idx<nsamp; idx++)
        adu_test[idx] = (double) data[(i*ncol+(j-1))*nsamp+idx];
      // look at pixel one column to the left
      corr = gsl_stats_correlation(adu_ref, 1, adu_test, 1, nsamp);
      printf("%g\n",corr);
      for (size_t idx=0; idx<nsamp; idx++)
        adu_test[idx] = (double) data[(i*ncol+(j+1))*nsamp+idx];
      // look at pixel one column to the right
      corr = gsl_stats_correlation(adu_ref, 1, adu_test, 1, nsamp);    
      printf("%g\n",corr);
    }
  }
  free(data);
  H5Fclose (file_id);
  return 0;
}
