#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#define baseline 5
struct data {
  gsl_matrix * stab_pxl; // matrix of stabilized observations from ROI
  double variance; // the CCD variance
};
int
residual_fct(const gsl_vector *log_par, void *data, 
	     gsl_vector * residual)
{
  gsl_matrix *obs = ((struct data *)data)->stab_pxl;
  double variance = ((struct data *)data)->variance;
  size_t n_pxl = obs->size1;
  size_t n_obs = obs->size2;
  size_t log_par_length = log_par->size;
  double auto_fluo = exp(gsl_vector_get(log_par,0));
  gsl_vector_const_view phi = gsl_vector_const_subvector(log_par, 1, n_pxl);
  gsl_vector_const_view f = gsl_vector_const_subvector(log_par, 1+n_pxl,n_obs-baseline);
  for (size_t pxl_idx=0; pxl_idx<n_pxl; pxl_idx++) {
    gsl_vector_const_view row = gsl_matrix_const_row(obs,pxl_idx);
    gsl_vector_view resid = gsl_vector_subvector(residual, pxl_idx*n_obs, n_obs);
    double phi_val = exp(gsl_vector_get(&phi.vector,pxl_idx));
    for (size_t k=0; k<baseline; k++) {
      double z = gsl_vector_get(&row.vector,k);
      double sF = 2*sqrt(auto_fluo+phi_val+variance);
      gsl_vector_set(&resid.vector,k,z-sF);
    }
    for (size_t k=baseline; k<n_obs; k++) {
      double z = gsl_vector_get(&row.vector,k);
      double s = exp(gsl_vector_get(&f.vector,k-baseline));
      double sF = 2*sqrt(auto_fluo+phi_val*s+variance);
      gsl_vector_set(&resid.vector,k,z-sF);
    }
  }
  return GSL_SUCCESS;
}
void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  
  fprintf(stderr, "iter %2zu: RSS = %.4f\n",
          iter,
          gsl_pow_2(gsl_blas_dnrm2(f)));
}
int main(int argc, char *argv[])
{
  static char usage[] = \
    "usage: %s -f --file=string -d --dset=string -g --CCD_gain=real ...\n"
    "          ... -v --CCD_variance=real -p --pixels=string \n\n"
    "  -f --file <character string>: data file name (e.g. Data_POMC.hdf5)\n"
    "  -d --dset <character string>: data set name (e.g. 'stack') shoud be a 3D object\n"
    "  -g --CCD_gain <positive real>: CCD gain (e.g. 0.14)\n"
    "  -v --CCV_variance <positive real>: CCD read-out variance (e.g. 290)\n"
    "  -p --pixels <character string>: name of a file containing the row and column\n"
    "       indexes of the pixels belonging to the ROI on two columns\n\n"
    " The program opens 'file' and gets dataset 'dset' (that should be a 3D array) containing\n"
    " ADUs. For each pixel i,j in the list specified by pixels at each time index k it\n"
    " stabilizes the variance:\n"
    "   adu_i,j,k -> z_i,j,k = 2*sqrt(adu_i,j,k/CCD_gain+CCD_variance)\n"
    " the variance of the z_i,j,k should then be UNIFORMLY 1. It then fits by nonlinear\n"
    " least-squares a model where the piecewise constant time course of the 'signal'\n"
    " of each pixel of the ROI is assumed to be the same. This 'signal' is then\n"
    " multiplied by a pixel specific factor (giving a predicted calcium induced)\n"
    " fluorescence and a common autofluorescence is added to each pixel.\n\n"
    " The program prints to the stderr fit diagnostics. It also prints to the stdout\n"
    " the 'signal' values together with the lower and upper limits of a 95% pointwise\n"
    " confidence intervals. Two blanck lines\n"
    " follow and, one as many columns as there are pixels in the ROI, the model\n"
    " predictions for each pixels are printed (a header gives the row and column\n"
    " index of each pixel). After two more blank lines the stabilized observations\n"
    " are printed in the same format.\n\n";
  char *filename;
  char dset[25];
  dset[0] = '/';
  dset[1] = '\0';
  double gain,variance;
  size_t *pxl_row,*pxl_col;
  size_t nb_pxl=0;
  {int opt;
    static struct option long_options[] = {
      {"file",required_argument,NULL,'f'},
      {"dset",required_argument,NULL,'d'},
      {"CCD_gain",required_argument,NULL,'g'},
      {"CCD_variance",required_argument,NULL,'v'},
      {"pixels",optional_argument,NULL,'p'},
      {"help",no_argument,NULL,'h'},
      {NULL,0,NULL,0}
    };
    char *dset_name;
    int long_index =0;
    while ((opt = getopt_long(argc,argv,
                              "hf:d:g:v:p:",
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
      case 'p':
      {
        char *fname = optarg;
        char line[BUFSIZ];
        FILE *fp;
        if (NULL == (fp = fopen(fname,"r"))) {
  	fprintf(stderr,"Could not open file %s\n",fname);
  	return -1;
        }
        while (fgets(line,BUFSIZ,fp))
  	nb_pxl++;
        rewind(fp);
        pxl_row = malloc(nb_pxl*sizeof(size_t));
        pxl_col = malloc(nb_pxl*sizeof(size_t));
        int row_idx,col_idx;
        size_t pxl_idx=0;
        while (fscanf(fp,"%d %d",&row_idx,&col_idx)==2) {
  	pxl_row[pxl_idx] = (size_t) row_idx;
  	pxl_col[pxl_idx] = (size_t) col_idx;
  	pxl_idx++;
        }
        fclose(fp);
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
  fprintf(stderr,"Fitting now piecewise constant signal within ROI for\n"
  	"dataset %s from file %s\n",dset,filename);
  fprintf(stderr,"The row and column indexes of the pixels making the\n"
  	"ROI are:\n");
  fprintf(stderr,"  (%d,%d)", (int) pxl_row[0],
  	  (int) pxl_col[0]);
  for (size_t pxl_idx=1; pxl_idx<nb_pxl; pxl_idx++)
    if (pxl_row[pxl_idx] > pxl_row[pxl_idx-1])
      fprintf(stderr,"\n  (%d,%d)", (int) pxl_row[pxl_idx],
  	    (int) pxl_col[pxl_idx]);
    else
      fprintf(stderr," (%d,%d)", (int) pxl_row[pxl_idx],
  	    (int) pxl_col[pxl_idx]);
  fprintf(stderr,"\nThe CCD gain and read-out variance are:\n"
  	"  %g and %g\n\n",gain,variance);
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
  gsl_matrix *pxl_data = gsl_matrix_alloc(nb_pxl,nsamp);
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++) {
    size_t i = pxl_row[pxl_idx];
    size_t j = pxl_col[pxl_idx];
    for (size_t k=0; k<nsamp; k++) {
      double adu = (double) data[(i*ncol+j)*nsamp+k];
      gsl_matrix_set(pxl_data,pxl_idx,k,adu);
    }
  }
  free(data);
  H5Fclose (file_id);
  gsl_vector *f = gsl_vector_calloc(nsamp);
  gsl_vector *phi = gsl_vector_calloc(nb_pxl);
  double b=100;
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++) {
    gsl_vector_const_view row = gsl_matrix_const_row(pxl_data,pxl_idx);
    gsl_vector_add(f,&row.vector);
  }
  gsl_vector_scale(f, 1.0/(double) nb_pxl);
  double baseline_mean=0;
  for (size_t k=0; k<baseline; k++)
    baseline_mean += gsl_vector_get(f,k);
  baseline_mean /= (double) baseline;
  gsl_vector_scale(f, 1.0/baseline_mean);
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++) {
    double sum=0;
    for (size_t k=0; k<baseline; k++)
      sum += gsl_matrix_get(pxl_data,pxl_idx,k);
    sum /= (double) baseline;
    gsl_vector_set(phi,pxl_idx,sum-b);
  }
  size_t n_par = 1+nb_pxl+(nsamp-baseline);
  gsl_vector *par = gsl_vector_alloc(n_par);
  gsl_vector_set(par,0,log(b));
  gsl_vector_view phiphi = gsl_vector_subvector(par,1,nb_pxl);
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++)
    gsl_vector_set(&phiphi.vector,pxl_idx,log(gsl_vector_get(phi,pxl_idx)));
  gsl_vector_view ff = gsl_vector_subvector(par,1+nb_pxl,nsamp-baseline);
  for (size_t k=baseline; k<nsamp; k++)
    gsl_vector_set(&ff.vector,k-baseline,log(gsl_vector_get(f,k)));
  // Stabilize variance of observations
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++) {
    for (size_t k=0; k<nsamp; k++) {
      // stabilize variance
      double z = 2*sqrt(gsl_matrix_get(pxl_data,pxl_idx,k)/gain+variance);
      gsl_matrix_set(pxl_data,pxl_idx,k,z);
    }
  }
  struct data d = { .stab_pxl = pxl_data,
                    .variance = variance};
  gsl_multifit_nlinear_fdf fdf;
  fdf.f = residual_fct;
  fdf.df = NULL;   /* set to NULL for finite-difference Jacobian */
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = nsamp*nb_pxl;
  fdf.p = n_par;
  fdf.params = &d;
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  /* allocate workspace with default parameters */
  gsl_multifit_nlinear_workspace *w = gsl_multifit_nlinear_alloc(T,
  							       &fdf_params,
  							       nsamp*nb_pxl,
  							       n_par);
  /* initialize solver with starting point */
  gsl_multifit_nlinear_init (par, &fdf, w);
  gsl_vector *residual = gsl_multifit_nlinear_residual(w);
  double chisq0;
  gsl_blas_ddot(residual, residual, &chisq0);
  const double xtol = 1e-10;
  const double gtol = 1e-10;
  const double ftol = 0.0;
  /* solve the system with a maximum of 200 iterations */
  int info;
  int status = gsl_multifit_nlinear_driver(200, xtol, gtol, ftol,
  					 callback, NULL, &info,
  					 w);
  
    /* compute covariance of best fit parameters */
  gsl_matrix *J = gsl_multifit_nlinear_jac(w);
  gsl_matrix *covar = gsl_matrix_alloc (1+nb_pxl+nsamp-baseline,
  				      1+nb_pxl+nsamp-baseline);
  gsl_multifit_nlinear_covar (J, 0.0, covar);
  
  /* compute final cost */
  double chisq;
  gsl_blas_ddot(residual, residual, &chisq);
  fprintf(stderr, "\n\nsummary from method '%s/%s'\n",
  	gsl_multifit_nlinear_name(w),
  	gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",
  	gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
  	(info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial RSS = %f\n", chisq0);
  fprintf(stderr, "final   RSS = %f\n", chisq);
  double dof = nsamp*nb_pxl - n_par;
  fprintf(stderr, "RSS/dof = %g\n", chisq / dof);
  fprintf (stderr, "status = %s\n\n", gsl_strerror (status));
  #define FIT(i) gsl_vector_get(w->x, i)
  #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
  printf("# Common background estimate: %g [%g,%g]\n\n",
         exp(FIT(0)),exp(FIT(0)-1.96*ERR(0)),
         exp(FIT(0)+1.96*ERR(0)));
  // Print f estimation as a first data set
  for (size_t k=0; k<baseline; k++)
    printf("%g %g %g\n",1.0,1.0,1.0);
  for (size_t k=baseline; k<nsamp; k++) {
    double est = FIT(1+nb_pxl+k-baseline);
    double err = ERR(1+nb_pxl+k-baseline);
    printf("%g %g %g\n",
  	 exp(est), exp(est-1.96*err),
  	 exp(est+1.96*err));
  }
  // Print intensity estimations one column per pixel
  printf("\n\n# Predicted variance stabilized signals");
  double auto_fluo = exp(FIT(0));
  printf("#");
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++)
    printf(" (%d,%d)",(int) pxl_row[pxl_idx],(int) pxl_col[pxl_idx]);
  printf("\n");
  for (size_t k=0; k<baseline; k++) {
    for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++)
      printf("%g ", 2*sqrt(exp(FIT(pxl_idx+1))+auto_fluo+variance));
    printf("\n");
  }
  for (size_t k=baseline; k<nsamp-baseline; k++) {
    for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++)
      printf("%g ", 2*sqrt(exp(FIT(pxl_idx+1))*exp(FIT(1+nb_pxl+k-baseline))+auto_fluo+variance));
    printf("\n");
  }
  // Print observations one column per pixel
  printf("\n\n# Observed stabilized signals");
  printf("#");
  for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++)
    printf(" (%d,%d)",(int) pxl_row[pxl_idx],(int) pxl_col[pxl_idx]);
  printf("\n");
  for (size_t k=0; k<nsamp; k++) {
    for (size_t pxl_idx=0; pxl_idx<nb_pxl; pxl_idx++)
      printf("%g ", gsl_matrix_get(pxl_data,pxl_idx,k));
    printf("\n");
  }
  free(pxl_row);
  free(pxl_col);
  gsl_matrix_free(pxl_data);
  gsl_vector_free(f);
  gsl_vector_free(phi);
  gsl_vector_free(par);
  gsl_multifit_nlinear_free (w);
  gsl_matrix_free(covar);  
  return 0;
}
