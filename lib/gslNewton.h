//
// Created by james on 10/11/16.
//

#ifndef CLIONPROJECTS_GSLNEWTON_H
#define CLIONPROJECTS_GSLNEWTON_H

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <stdlib.h>

// ----------- Global Variables ------- //

int param_len_newt ;

// --------------------------------- //

struct func_Params_newt{
    double *d;
};

double ratFun2_newton (double x, void *params)
{
  /*
   * ---------------------------------------------
   * Rational Function of the form :
   * f(x) = a1 / ( x + b1 ) + a2 / (x + b2 )
   * ---------------------------------------------
   */
  struct func_Params_newt *p = (struct func_Params_newt *) params;
  double sum = 0;
  int i, j = 0;
  for(i = 0 ; i < param_len_newt/2 ; i++){
      sum = sum + (p->d[j] / (x + p->d[j+1]));
      j = j + 2;
  }
  return sum;
}

double ratFun2_deriv (double x, void *params)
{
  /*
   * ---------------------------------------------
   * Rational Function of the form :
   * f'(x) = -1 * a1 / ( x + b1 )^2 + (-1)*a2 / (x + b2 )^2
   * ---------------------------------------------
   */
  struct func_Params_newt *p = (struct func_Params_newt *) params;
  double sum = 0;
  int i, j = 0;
  for(i = 0 ; i < param_len_newt/2 ; i++){
      sum = sum + ((-1)*p->d[j] / ((x + p->d[j+1])*(x + p->d[j+1])));
      j = j + 2;
  }
  return sum;
}

void ratFun2_fdf (double x, void *params,
               double *y, double *dy)
{

  struct func_Params_newt *p = (struct func_Params_newt *) params;

  // --------------- y --------- //

  double sum_y = 0;
  int i, j = 0;
  for(i = 0 ; i < param_len_newt/2 ; i++){
      sum_y = sum_y + (p->d[j] / (x + p->d[j+1]));
      j = j + 2;
  }

  // ------------- dy ----------- //
  double sum_dy = 0;
  j = 0;
  for(i = 0 ; i < param_len_newt/2 ; i++){
      sum_dy = sum_dy + ((-1)*p->d[j] / ((x + p->d[j+1])*(x + p->d[j+1])));
      j = j + 2;
  }

  *y = sum_y;
  *dy = sum_dy;
}

long double getRoot2_newton(double *pdata, int n, double x)
{
    int status;
    param_len_newt = n;
    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0;
    gsl_function_fdf FDF;
    struct func_Params_newt params;
    params.d = pdata;

    FDF.f = &ratFun2_newton;
    FDF.df = &ratFun2_deriv;
    FDF.fdf = &ratFun2_fdf;
    FDF.params = &params;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_root_fdfsolver_set (s, &FDF, x);
    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-3);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fdfsolver_free (s);
    // printf("From C : %f\n", x);
    return x;
}

double ratFun_newton (double x, void *params)
{
  /*
   * ---------------------------------------------
   * Rational Function of the form :
   * f(x) = a1 / ( x + b1 ) + a2 / (x + b2 ) + c
   * ---------------------------------------------
   */
   int tmpLen;
   tmpLen = param_len_newt - 1;
   struct func_Params_newt *p = (struct func_Params_newt *) params;
   double sum = 0;
   int i, j = 0;
   for(i = 0 ; i < tmpLen/2 ; i++){
       sum = sum + (p->d[j] / (x + p->d[j+1]));
       j = j + 2;
   }
   sum = sum + p->d[tmpLen];
   return sum;
}

double ratFun_deriv (double x, void *params)
{
  /*
   * ------------------------------------------------------
   * Rational Function of the form :
   * f'(x) = -1 * a1 / ( x + b1 )^2 + (-1)*a2 / (x + b2 )^2
   * ------------------------------------------------------
   */
   int tmpLen;
   tmpLen = param_len_newt - 1;
   struct func_Params_newt *p = (struct func_Params_newt *) params;
   double sum = 0;
   int i, j = 0;
   for(i = 0 ; i < tmpLen/2 ; i++){
       sum = sum + ((-1)*p->d[j] / ((x + p->d[j+1])*(x + p->d[j+1])));
       j = j + 2;
   }
   return sum;
}

void ratFun_fdf (double x, void *params,
               double *y, double *dy)
{

  // --------------- y --------- //
  double sum_y = 0;
  int tmpLen;
  tmpLen = param_len_newt - 1;
  struct func_Params_newt *p = (struct func_Params_newt *) params;
  double sum = 0;
  int i, j = 0;
  for(i = 0 ; i < tmpLen/2 ; i++){
      sum_y = sum_y + (p->d[j] / (x + p->d[j+1]));
      j = j + 2;
  }
  sum_y = sum_y + p->d[tmpLen];
  *y = sum_y;

  // ------------- dy ----------- //

  tmpLen = param_len_newt - 1;
  double sum_dy = 0;
  i, j = 0;
  for(i = 0 ; i < tmpLen/2 ; i++){
      sum_dy = sum_dy + ((-1)*p->d[j] / ((x + p->d[j+1])*(x + p->d[j+1])));
      j = j + 2;
  }
  *dy = sum_dy;
}

long double getRoot_newton(double *pdata, int n, double x)
{
    int status;
    param_len_newt = n;
    int iter = 0, max_iter = 100;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double x0;
    gsl_function_fdf FDF;
    struct func_Params_newt params;
    params.d = pdata;

    FDF.f = &ratFun_newton;
    FDF.df = &ratFun_deriv;
    FDF.fdf = &ratFun_fdf;
    FDF.params = &params;

    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc(T);
    gsl_root_fdfsolver_set (s, &FDF, x);
    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 1e-3);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fdfsolver_free (s);
    // printf("From C : %f\n", x);
    return x;
}
#endif //CLIONPROJECTS_GSLNEWTON_H
