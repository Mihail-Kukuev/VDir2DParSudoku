#ifndef _TASK_DEF_H

#define _TASK_DEF_H


#include <stdio.h>
#include <math.h>


#define L1 1
#define L2 1
#define T  1


// Test solution
double u ( double x1, double x2, double t );


// Conditions
double u0 ( double x1, double x2 ); // t = 0

double m1 ( double x2, double t );  // x1 = 0
double m2 ( double x2, double t );  // x1 = L1
double m3 ( double x1, double t );  // x2 = 0
double m4 ( double x1, double t );  // x2 = L2


// Coefficients
double k1 ( double x1, double x2, double t );
double k2 ( double x1, double x2, double t );

double f ( double x1, double x2, double t );

#endif
