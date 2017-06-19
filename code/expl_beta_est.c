#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#define PRINT(what) {							\
  fprintf(stderr, "summary from method '%s/%s'\n",			\
	  gsl_multifit_nlinear_name(w_##what),				\
          gsl_multifit_nlinear_trs_name(w_##what));			\
  fprintf(stderr, "number of iterations: %zu\n",			\
          gsl_multifit_nlinear_niter(w_##what));			\
  fprintf(stderr, "function evaluations: %zu\n", fdf_##what.nevalf);	\
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf_##what.nevaldf);	\
  fprintf(stderr, "reason for stopping: %s\n",				\
          (info_##what == 1) ? "small step size" : "small gradient");	\
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0_##what));	\
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq_##what));		\
  fprintf (stderr, "beta = %.5f +/- %.5f\n",				\
	   gsl_vector_get(w_##what->x, 0),				\
	   1.96*sqrt(gsl_matrix_get(covar_##what,0,0)));		\
  fprintf (stderr, "status = %s\n", gsl_strerror (status_##what));	\
}
struct data {
  size_t n; // number of observations
  double * t; // pointer to an array (of size n) holding t_1 and t_2
  double * y; // pointer to an array (of size n) holding y_1 and y_2
  double F_inf,Delta,t_onset; // the known model parameters
};
int
tilde_f (const gsl_vector * x, void *data, 
	 gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *t = ((struct data *)data)->t;
  double F_inf = ((struct data *)data)->F_inf;
  double Delta = ((struct data *)data)->Delta;
  double t_onset = ((struct data *)data)->t_onset;
  double beta = gsl_vector_get (x, 0);

  for (size_t i = 0; i < n; i++) {
    /* Model Yi = F_inf + H(t-t_onset) * Delta * exp(-beta * t) */
    double Yi = F_inf;
    double delta_t = t[i]-t_onset;
    if (delta_t >= 0.) 
      Yi += Delta * exp (-beta * delta_t);
    gsl_vector_set (f, i, Yi - y[i]);
  }

  return GSL_SUCCESS;
}
int
tilde_df (const gsl_vector * x, void *data, 
	  gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *t = ((struct data *)data)->t;
  double F_inf = ((struct data *)data)->F_inf;
  double Delta = ((struct data *)data)->Delta;
  double t_onset = ((struct data *)data)->t_onset;
  double beta = gsl_vector_get (x, 0);

  for (size_t i = 0; i < n; i++) {
    /* Jacobian matrix J(i,j) = dfi / dxj, */
    /* where fi = (Yi - yi)/sigma[i],      */
    /*       Yi = F_inf + H(t-t_onset) * Delta * exp(-beta * t)  */
    /* and xj is the model parameter beta */
    double delta_t = t[i] - t_onset;
    if (delta_t >= 0) {
      double e = exp(-beta * delta_t);
      gsl_matrix_set (J, i, 0, -delta_t * Delta * e);
    } else {
      gsl_matrix_set (J, i, 0, 0);
    }
  }
  return GSL_SUCCESS;
}
int
hat_f (const gsl_vector * x, void *data, 
       gsl_vector * f)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *t = ((struct data *)data)->t;
  double F_inf = ((struct data *)data)->F_inf;
  double Delta = ((struct data *)data)->Delta;
  double t_onset = ((struct data *)data)->t_onset;
  double beta = gsl_vector_get (x, 0);

  for (size_t i = 0; i < n; i++) {
    /* Model Yi = sqrt(F_inf + H(t-t_onset) * Delta * exp(-beta * t)) */
    double Yi = F_inf;
    double delta_t = t[i]-t_onset;
    if (delta_t >= 0) 
      Yi += Delta * exp (-beta * delta_t);
    gsl_vector_set (f, i, sqrt(Yi) - sqrt(y[i]));
  }

  return GSL_SUCCESS;
}
int
hat_df (const gsl_vector * x, void *data, 
	gsl_matrix * J)
{
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *t = ((struct data *)data)->t;
  double F_inf = ((struct data *)data)->F_inf;
  double Delta = ((struct data *)data)->Delta;
  double t_onset = ((struct data *)data)->t_onset;
  double beta = gsl_vector_get (x, 0);

  for (size_t i = 0; i < n; i++) {
    /* Jacobian matrix J(i,j) = dfi / dxj, */
    /* where fi = (Yi - yi)/sigma[i],      */
    /*       Yi = F_inf + H(t-t_onset) * Delta * exp(-beta * t)  */
    /* and xj is the model parameter beta */
    double delta_t = t[i] - t_onset;
    if (delta_t >= 0) {
      double e = exp(-beta * delta_t);
      double Yi = F_inf+Delta*e;
	gsl_matrix_set (J, i, 0, -0.5 * delta_t * Delta * e / sqrt(Yi));
    } else {
      gsl_matrix_set (J, i, 0, 0);
    }
  }
  return GSL_SUCCESS;
}
void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  
fprintf(stderr, "iter %2zu: beta = %.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_blas_dnrm2(f));
}
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
  double y[n_obs], t[n_obs];
  for (size_t i=0; i<n_obs; i++) {
    t[i] = t_0 + i*delta_t;
    double ideal = F_inf;
    if (t[i] >= t_onset) ideal += Delta*exp(-beta*(t[i]-t_onset));
    y[i] = gsl_ran_poisson(r, ideal);
  }
  struct data d = { .n = n_obs, .t = t, .y = y,
  		  .F_inf = F_inf, .Delta = Delta,
  		  .t_onset = t_onset};
  double x_init[1] = { beta }; /* starting values */
  gsl_vector_view x = gsl_vector_view_array (x_init, 1);
  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  gsl_multifit_nlinear_fdf fdf_tilde;
  fdf_tilde.f = tilde_f;
  fdf_tilde.df = tilde_df;   /* set to NULL for finite-difference Jacobian */
  fdf_tilde.fvv = NULL;     /* not using geodesic acceleration */
  fdf_tilde.n = n_obs;
  fdf_tilde.p = 1;
  fdf_tilde.params = &d;
  const gsl_multifit_nlinear_type *T_tilde = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params_tilde =
    gsl_multifit_nlinear_default_parameters();
  /* allocate workspace with default parameters */
  gsl_multifit_nlinear_workspace *w_tilde = gsl_multifit_nlinear_alloc (T_tilde,
  								      &fdf_params_tilde,
  								      n_obs, 1);
  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_init (&x.vector, &fdf_tilde, w_tilde);
  gsl_vector *f_tilde = gsl_multifit_nlinear_residual(w_tilde);
  double chisq0_tilde;
  gsl_blas_ddot(f_tilde, f_tilde, &chisq0_tilde);
  /* solve the system with a maximum of 20 iterations */
  int info_tilde;
  int status_tilde = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol,
  					       callback, NULL, &info_tilde,
  					       w_tilde);
  
    /* compute covariance of best fit parameters */
  gsl_matrix *J_tilde = gsl_multifit_nlinear_jac(w_tilde);
  gsl_matrix *covar_tilde = gsl_matrix_alloc (1, 1);
  gsl_multifit_nlinear_covar (J_tilde, 0.0, covar_tilde);
  
  /* compute final cost */
  double chisq_tilde;
  gsl_blas_ddot(f_tilde, f_tilde, &chisq_tilde);
  PRINT(tilde)

  gsl_multifit_nlinear_fdf fdf_hat;
  fdf_hat.f = hat_f;
  fdf_hat.df = hat_df;   /* set to NULL for finite-difference Jacobian */
  fdf_hat.fvv = NULL;     /* not using geodesic acceleration */
  fdf_hat.n = n_obs;
  fdf_hat.p = 1;
  fdf_hat.params = &d;
  const gsl_multifit_nlinear_type *T_hat = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params_hat =
    gsl_multifit_nlinear_default_parameters();
  /* allocate workspace with default parameters */
  gsl_multifit_nlinear_workspace *w_hat = gsl_multifit_nlinear_alloc (T_hat,
  								    &fdf_params_hat,
  								    n_obs, 1);
  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_init (&x.vector, &fdf_hat, w_hat);
  gsl_vector *f_hat = gsl_multifit_nlinear_residual(w_hat);
  double chisq0_hat;
  gsl_blas_ddot(f_hat, f_hat, &chisq0_hat);
  /* solve the system with a maximum of 20 iterations */
  int info_hat;
  int status_hat = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol,
  					       callback, NULL, &info_hat,
  					       w_hat);
  
    /* compute covariance of best fit parameters */
  gsl_matrix *J_hat = gsl_multifit_nlinear_jac(w_hat);
  gsl_matrix *covar_hat = gsl_matrix_alloc (1, 1);
  gsl_multifit_nlinear_covar (J_hat, 0.0, covar_hat);
  
  /* compute final cost */
  double chisq_hat;
  gsl_blas_ddot(f_hat, f_hat, &chisq_hat);
  PRINT(hat)

  gsl_rng_free(r);
  gsl_multifit_nlinear_free (w_tilde);
  gsl_multifit_nlinear_free (w_hat);
  gsl_matrix_free (covar_tilde);
  gsl_matrix_free (covar_hat);
  return 0;
}
