#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s -f --file=string -d --dset=string -g --CCD_gain=real ...\n"
    "          ... -v --CCD_variance=real [-t --threshold=real] [-m --return_matrix]\n\n"
    "  -f --file <character string>: data file name (e.g. Data_POMC.hdf5)\n"
    "  -d --dset <character string>: data set name (e.g. 'stack') shoud be a 3D object\n"
    "  -g --CCD_gain <positive real>: CCD gain (e.g. 0.14)\n"
    "  -v --CCV_variance <positive real>: CCD read-out variance (e.g. 290)\n"
    "  -m --return_matrix if set a matrix whose values at pixel i,j are log(1-F(rss_i,j))\n"
    "     where F is the CDF of a chi square distribution with nsamp-1 degrees of freedome\n"
    "     nsamp being the number of measurements (exposures) made on each pixel. If unset\n"
    "     (default) the list of pixel in the ROI is returned as (i1,j1),(i2,j2),... \n"
    "  -t --threshold <negative real>: real number (e.g. -300 default), the level\n"
    "     under which the  log(1-F(rss_i,j)) is considered significant\n\n"
    " The program opens 'file' and gets dataset 'dset' (that should be a 3D array) containing\n"
    " ADUs. For each pixel i,j at each time index k it stabilizes the variance:\n"
    "   adu_i,j,k -> z_i,j,k = 2*sqrt(adu_i,j,k/CCD_gain+CCD_variance)\n"
    " the variance of the z_i,j,k should then be UNIFORMLY 1. It then computes the mean value\n"
    " along time for each pixel:\n"
    "   z_i,j_m = (z_i,j,1 + ... + z_i,j,K)/K\n"
    " and the residual sum of squares statistic:\n"
    "   rss_i,j = ((z_i,j,1-z_i,j_m)^2+...+(z_i,j,K-z_i,j_m)^2))\n"
    " under the null hypothesis (no signal in the pixel) this statistic should have a chi-square\n"
    " distribution with K-1 degrees of freedom. A matrix whose element i,j is:\n"
    "   log(1-F(rss_i,j))\n"
    " where F is the CDF of a chi-square distribution with K-1 degrees of freedom is constructed.\n"
    " If parameter '--return_matrix' is set (not the default) this matrix is printed to the stdout,\n"
    " otherwise the value of the F(rss_i,j) is compared to the Bonferroni corrected threshold set\n"
    " by '--threshold' (default 0.99) and if the value is larger, 'i j' is printed to the stdout.\n"
    " This is done for each pixel and the printed pairs are separated by '\n'.\n\n";
  char *filename;
  char dset[25];
  dset[0] = '/';
  dset[1] = '\0';
  double threshold = -300.;
  double gain,variance;
  int print_matrix=0;
  {int opt;
    static struct option long_options[] = {
      {"file",required_argument,NULL,'f'},
      {"dset",required_argument,NULL,'d'},
      {"CCD_gain",required_argument,NULL,'g'},
      {"CCD_variance",required_argument,NULL,'v'},
      {"threshold",optional_argument,NULL,'t'},
      {"return_matrix",no_argument,NULL,'m'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    char *dset_name;
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hmf:d:g:v:t:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'f':
      {
        filename = optarg;
      }
      break;
      case 'd':
      {
        dset_name = optarg;
      }
      break;
      case 'g':
      {
        gain = (double) atof(optarg);
        if (gain <= 0.0) {
          fprintf(stderr,"CCD gain must be > 0\n");
          return -1;
        }
      }
      break;
      case 'v':
      {
        variance = (double) atof(optarg);
        if (variance < 0.0) {
          fprintf(stderr,"CCD variance must be >= 0\n");
          return -1;
        }
      }
      break;
      case 't':
      {
        threshold = (double) atof(optarg);
        if (0.0 <= threshold) {
          fprintf(stderr,"threshold should be\n");
          return -1;
        }
      }
      break;
      case 'm':
      {
        print_matrix=1;
      }
      break;
      case 'h': printf(usage,argv[0]);
        return -1;
      default : fprintf(stderr,usage,argv[0]);
        return -1;
      }
    }
    strcat(dset,dset_name);
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
  gsl_matrix *pRSS = gsl_matrix_alloc(nrow,ncol);
  for (size_t i=0; i<nrow; i++) {
    for (size_t j=0; j<ncol; j++) {
      double z[nsamp];
      for (size_t k=0; k<nsamp; k++) {
        double adu = (double) data[(i*ncol+j)*nsamp+k];
        z[k] = 2.0*sqrt(adu/gain+variance);
      }
      double rss = gsl_stats_variance(z,1,nsamp)*((double)nsamp-1.0); 
      gsl_matrix_set(pRSS,i,j,log(gsl_cdf_chisq_Q(rss, (double)nsamp-1.0)));
    }
  }
  if (print_matrix) {
    printf("# Matrix with %d rows and %d columns\n",
  	 (int) nrow, (int) ncol);
    for (size_t i=0; i<nrow; i++) {
      printf("%g",gsl_matrix_get(pRSS,i,0));
      for (size_t j=1; j<ncol; j++)
        printf(", %g",gsl_matrix_get(pRSS,i,j));
      printf("\n");
    }
  } else {
    size_t n_pixel=0;
    for (size_t i=0; i<nrow; i++) {
      for (size_t j=1; j<ncol; j++) {
        double stat = gsl_matrix_get(pRSS,i,j);
        if (stat <= threshold) {
  	// Pixel part of ROI
  	if (n_pixel==0)
  	  printf("%d %d", (int) i, (int) j);
  	else
  	  printf("\n%d %d", (int) i, (int) j);
  	n_pixel++;
        }
      }
    }
    if (n_pixel > 0)
      printf("\n");
    fprintf(stderr,"%d pixels in the ROI\n", (int) n_pixel);
  }
  free(data);
  gsl_matrix_free(pRSS);
  H5Fclose (file_id);  
  return 0;
}
