
/*
    Linear least squares diagnostics
    Copyright (C) 2024  T.P. Leão

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    (For) a copy of the GNU General Public License ... see <https://www.gnu.org/licenses/>.
*/

/*** RUN THE CODE WITH THESE COMMANDS IN TERMINAL
gcc -Wall -I/usr/local/include -c linreg.c
gcc -L/usr/local/lib linreg.o -lgsl -lgslcblas -lm -o linreg

LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH
./linreg
**/

// Code starts
#include <stdio.h> 		//Standard input-output functions
#include <string.h>		//String functions
#include <gsl/gsl_fit.h> 	//GSL linear least-squares functions 
#include <gsl/gsl_statistics.h> //GSL statistics functions
#include <gsl/gsl_randist.h> 	//GSL distributions functions
#include <gsl/gsl_cdf.h> 	//GSL cumulative distributions functions


// Function for significance level of parameters
const char * significance(double a)
{
 static char s[6];
  	if (a > 0.05){strcpy(s, "ns"); } 
  	else if (a <= 0.05 && a > 0.01) {strcpy(s, "*");  } 
 	else if (a <= 0.01 && a > 0.001) {strcpy(s, "**");} 
 	else if (a <=  0.001) {strcpy(s, "***");} 
 	else  {strcpy(s, "error");} 
	return s;
}


// Main program

int
main (void)
{
  int i, n = 15;
  
 // Data - add your dataset here
  double x[15] = { 20, 16, 20, 18, 17, 16, 15, 17, 15, 16, 15, 17, 16, 17, 14 };
  double y[15] = { 89 ,72 ,93 ,84, 81, 75, 70, 82, 69, 83, 80, 83, 81, 84, 76 };
  
  
  double yp[n], res[n];
  double c0, c1, cov00, cov01, cov11, sumsq, cov, r, spear, se_c0, se_c1, t_c0, t_c1, p_c0, p_c1;


//Functions

//Linear regression
gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);                   

//Covariance 
cov = gsl_stats_covariance(x, 1, y, 1, n);

//Correlation
r = gsl_stats_correlation(x, 1, y, 1, n);

//Standard error of parameter estimates
se_c0 = sqrt(cov00);
se_c1 = sqrt(cov11);

//t-value
t_c0 = c0/se_c0;
t_c1 = c1/se_c1;

//p-value
p_c0 = 2*gsl_cdf_tdist_P(-fabs(t_c0), n-2);
p_c1 = 2*gsl_cdf_tdist_P(-fabs(t_c1), n-2);

//Spearman correlation
double work[2*n];
spear = gsl_stats_spearman(x, 1, y, 1, n, work);

//Regression diagnostics

  printf ("###################### REGRESSION DIAGNOSTICS ########################## \n");
   
  printf ("# Model: y = c0 + c1 x \n");
  printf ("# Best fit: y = %g + %g x\n", c0, c1);
  printf ("# Sum of squares of residuals:  %g \n", sumsq);
  printf ("# Standard error of estimates: c0 = %g, c1 = %g  \n", se_c0, se_c1);
  printf ("# t-value of estimates: c0 = %g, c1 = %g  \n", t_c0, t_c1);
  printf ("# Associated probability values: c0 = %g, c1 = %g  \n", p_c0, p_c1);
  printf ("# Associated significance values: c0 = %s, " , significance(p_c0));
  printf (" c1 = %s    \n" , significance(p_c1));
  printf ("# Interpretation: ns = not significant, * = p <= 0.05, ** = p <= 0.01, *** = p <= 0.001 \n");
  
  printf ("# Covariance matrix:\n");
  printf ("# [ %g, %g\n#   %g, %g]\n", cov00, cov01, cov01, cov11);
  printf ("# Covariance COV:  %.2g \n", cov);
  printf ("# Correlation r: %.2g \n", r);
  printf ("# Coeficient of determination r2: %.2g \n", pow(r,2));
  printf ("# Spearman correlation:  %.2g \n", spear);


// Prints the dataset to terminal, incluind predicted values
 printf ("############################ DATASET ###################################### \n");
 

printf ("# Dataset: \n");
printf ("x  y  yp residuals \n");

  for (i = 0; i < n; i++){
   yp[i] = c0 + c1*x[i];	
   res[i] = yp[i] - y[i];
    printf ("%g %g %g %g \n",x[i], y[i], yp[i], res[i]);
}

  printf ("\n");

/* Use this code fragment to extrapolate, including calculated error bars
  for (i = -30; i < 130; i++)
    {
      double xf = x[0] + (i/100.0) * (x[n-1] - x[0]);
      double yf, yf_err;

      gsl_fit_linear_est (xf,
                          c0, c1,
                          cov00, cov01, cov11,
                          &yf, &yf_err);

      printf ("fit: %g %g\n", xf, yf);
      printf ("hi : %g %g\n", xf, yf + yf_err);
      printf ("lo : %g %g\n", xf, yf - yf_err);
    }
 */
   printf ("############################ END ###################################### \n");

  return 0;
}
