#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#define FILE "Data_POMC.hdf5"
#define DSET "stack"

int main(int argc, char *argv[])
{
  if (argc != 5) {
    fprintf(stderr,"Expecting four arguments\n");
    return -1;
  }
  size_t first_row = atoi(argv[1]);
  size_t last_row = atoi(argv[2]);
  size_t first_col = atoi(argv[3]);
  size_t last_col = atoi(argv[4]);
  // Open FILE
  hid_t file_id = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
  // Get dimensions of 3D object contained in DSET
  hsize_t dims[3];
  H5LTget_dataset_info(file_id,DSET,dims,NULL,NULL);
  size_t nrow = (size_t) dims[0];
  size_t ncol = (size_t) dims[1];
  size_t nsamp = (size_t) dims[2];
  // Read DSET in DATA
  int data[nrow][ncol][nsamp];
  H5LTread_dataset_int(file_id,"/stack",data);
  // Find out the smallest and largest observations
  // in the selected part of DSET
  int adu_max=0;
  int adu_min=10000;
  for (size_t i=first_row; i<last_row; i++) {
    for (size_t j=first_col; j<last_col; j++)  {
      for (size_t k=0; k<nsamp; k++) {
        double adu = data[i][j][k];
        if (adu < adu_min) adu_min=adu;
        if (adu > adu_max) adu_max=adu;
      }
    }
  }
  // Print some info to STDOUT
  printf("# Data set stack from file: %s\n",FILE);
  printf("# Data set dimensions: (%d,%d,%d)\n",nrow,ncol,nsamp);
  printf("# Using rows from %d (inclusive) to %d (exclusive)\n",first_row,last_row);
  printf("# Using columns from %d (inclusive) to %d (exclusive)\n",first_col,last_col);
  printf("# Minimal ADU in this range: %d; maximal value: %d\n",adu_min,adu_max);
  double adu_delta = (double) (adu_max-adu_min);
  // Write the DATA in a 2 columns format with time in the first and normalized
  // ADU in the second
  for (size_t i=first_row; i<last_row; i++) {
    double y_min = i-first_row;
    for (size_t j=first_col; j<last_col; j++)  {
      for (size_t k=0; k<nsamp; k++) {
        double adu = (data[i][j][k]-adu_min)/adu_delta+y_min;
        printf("%g %g\n",((double) k/nsamp+j-first_col),adu);
      }
      printf("\n");
    }
  }
  H5Fclose (file_id);
  return 0;

}
