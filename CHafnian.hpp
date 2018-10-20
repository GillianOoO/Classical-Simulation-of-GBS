//
//  CHafnian.hpp
//  GBSNoncollision
//
//  Created by WuBujiao on 8/15/18.
//  Copyright © 2018 WuBujiao. All rights reserved.
//

#ifndef CHafnian_hpp
#define CHafnian_hpp

#include <complex.h>  /* Change this <complex.h> into </path/complex.h>*/
#include <assert.h>
#include </usr/local/opt/openblas/include/lapacke.h>
#include </usr/local/include/omp.h>
//#include <lapacke.h>
//#include <omp.h>
#include <math.h>


typedef double complex telem;
typedef unsigned char sint;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void haf (telem mat[], int n, double res[]);
telem lhafnian (telem mat[], int n);// return hafnian of mat[n*n]
void evals (double complex z[], double complex vals[], int n);
telem hafnian_loops(telem *mat, int n);

#endif