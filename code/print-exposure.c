#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <getopt.h>
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s [--group=string] [--frame=idx] [--file=string]\n\n"\
    "  --group <character string>: the name of the group (eg, 10ms the default)\n"
    "  --frame <positive integer>: the frame index (from 0 to 99, default at 0)\n"
    "  --file <character string>: data file name (default CCD_calibration.hdf5)\n\n"
    " The program opens 'file', gets to dataset 'group/stack' and prints to\n"
    " the stdout the exposure whose index is given by frame.\n\n";
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
      {"frame",optional_argument,NULL,'i'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    int long_index =0;
    int file_name_set = 0, group_name_set = 0;
    while ((opt = getopt_long(argc,argv,
                              "hf:g:i:",
                              long_options,       
                              &long_index)) != -1) {
      switch(opt) {
      case 'i':
      {
        frame = (size_t) atoi(optarg);
        if (frame > 99)
        {
  	fprintf(stderr,"frame should between 0 and 99 inclusive.\n");
  	return -1;
        }  
      }
      break;
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
  for (size_t i=0; i<nrow; i++) {
    for (size_t j=0; j<ncol; j++)
      printf("%d ",data[(i*ncol+j)*nsamp+frame]);
    printf("\n");
  }
  free(data);
  H5Fclose (file_id); 
  return 0;
}
