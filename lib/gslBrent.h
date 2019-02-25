//
// Created by james on 8/23/16.
//

#ifndef EM_DCOP_GSLBRENT_H
#define EM_DCOP_GSLBRENT_H


#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <stdlib.h>


// ----------- Global Variables ------- //

int param_len ;

// --------------------------------- //

struct func_Params{
    double *d;
};

double ratFun(double x, void *params){

    /*
     * ---------------------------------------------
     * Rational Function of the form :
     * f(x) = a1 / ( x + b1 ) + a2 / (x + b2 ) + c
     * ---------------------------------------------
     */

    int tmpLen;
    tmpLen = param_len - 1;
    struct func_Params *p = (struct func_Params *) params;
    double sum = 0;
    int i, j = 0;
    for(i = 0 ; i < tmpLen/2 ; i++){
        sum = sum + (p->d[j] / (x + p->d[j+1]));
        j = j + 2;
    }
    sum = sum + p->d[tmpLen];
    return sum;
}


double ratFun2(double x, void *params){

    /*
     * ---------------------------------------------
     * Rational Function of the form :
     * f(x) = a1 / ( x + b1 ) + a2 / (x + b2 )
     * ---------------------------------------------
     */

    struct func_Params *p = (struct func_Params *) params;
    double sum = 0;
    int i, j = 0;
    for(i = 0 ; i < param_len/2 ; i++){
        sum = sum + (p->d[j] / (x + p->d[j+1]));
        j = j + 2;
    }
    return sum;
}

//(*nq)[n]
double getRoot (double *pdata, int n, double l, double h)
{

    int status;
    param_len = n;
    int iter = 0, max_iter = 500;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double x_lo, x_hi;
    x_lo = l;
    x_hi = h;
    double r = 0;
    gsl_function F;
    struct func_Params params;
    params.d = pdata;
    F.function = &ratFun;
    F.params = &params;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    //printf("From C : %f\n", r);
    return r;
}


long double getRoot2 (double *pdata, int n, double l, double h)
{

    int status;
    param_len = n;
    int iter = 0, max_iter = 500;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double x_lo, x_hi;
    x_lo = l;
    x_hi = h;
    double r = 0;
    gsl_function F;
    struct func_Params params;
    params.d = pdata;
    F.function = &ratFun2;
    F.params = &params;
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);
    // printf("From C : %f\n", r);
    return r;
}

#endif //EM_DCOP_GSLBRENT_H
