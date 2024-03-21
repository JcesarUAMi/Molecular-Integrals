#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <CL/sycl.hpp>
#include "BoysFunction.h"

using namespace std;
using namespace sycl;


/*template<typename U> U Boys<double>::boysUpward(U x, int n) {

  U result, value;
  value = 1.0/(2.0*x);
  result = SplinesBoys<double>(9, x);
  for (int i = 50; i < n; i++) {
    result *= value * (2*i + 1);
    result -= value * exp(-x);
  }
  return result;

}

template<typename T> T Boys<double>::boysDownward(T integral, T x, int m, int p) {

  T result;
  int n, j;
  j = m - p;
  result = integral;
  for (int i = 1; i <= j; i++) {
    n = m - i;
    result *= 2 * x;
    result += exp(-x);
    result /= (2 * n + 1);
  }

  return result;
}


template <typename T> T Boys::SplinesBoys(int li, T xj) {

  int ind_j, val;
  T fx, xi, sxi;

  ind_j = xj * 10000;
  xi = xj - ind_j * 0.0001;
	sxi = xi*xi;

	if (li > 2) {
		val = li/5;	
  	fx = btv[val].a[ind_j] + btv[val].b[ind_j]*xi + btv[val].c[ind_j]*sxi + btv[val].d[ind_j]*sxi*xi;
	} else 
  	fx = btvNum[li].a[ind_j] + btvNum[li].b[ind_j]*xi + btvNum[li].c[ind_j]*sxi + btvNum[li].d[ind_j]*sxi*xi;
	
  return fx;
}

template <typename T> T Boys::boysSeriesInferior (int n, T x) {
  T numerator, denominator, result;
  numerator = tgamma(n+0.5);
  result = 0.0;
  for (int i=0; i<50; i++) {
    denominator = tgamma(n+i+1.5);
    result += pow(x, i)*numerator/denominator;
  }
  result *= 0.5 * exp(-x);

  return result;
}

template <typename T> T Boys::boysZeroDegree(T x) {
  
  T result, sqX;

  if (x == 0.0)
   result = 1.0;
  else { 
    sqX = sqrt(x); 
    if (sqX >= 6.0) 
      result = sqPi/sqX;
    else
      result = erf(sqX)*sqPi/sqX;
  }

  return result;
}

template <typename T> T Boys::factorialDouble(int n) {
  int value;
  int result;
  result = 1;
  if (n > 0) {
    value = n/2 - 1;
    for (int i=0; i<=value; i++)
      result *= (n - 2*i);
  }

  return result;
}

template <typename T> T Boys::BoysFunction(int n, T xj) {

  int k, ni, i, j; 
  div_t m; 
  T integral, fx;


	if (xj < 1E-12)
		xj = 0.0;		

	if (xj == 0) 
  	fx = double(1.0 / (2.0*n + 1.0));
 	else if (xj < 1E-4) {
    ni = n/5;
    integral = SplinesBoys<double>(ni,xj);
    ni = (ni + 1) * 5;
    fx = boysDownward<double>(integral, xj, ni, n);
  }	else if (xj > 30.0) {
    fx = factorialDouble<double>(2*n-1)/(pow(2.0,n+1));
    fx *= sqrt(M_PI/pow(xj,2*n+1));
  } else if (n < 4) 
		fx = SplinesBoys<double>(n-1, xj);
	else if (n > 3 && n < 50) { 
    m = div(n, 5);
		i = m.rem;
		j = m.quot;
    if (i == 0) 
      fx = SplinesBoys<double>(n-1, xj);
    else {
      integral = SplinesBoys<double>(n, xj);
      k = (j + 1) * 5;
      fx = boysDownward<double>(integral, xj, k, n);
		}  
  } else 
		fx = boysUpward<double>(xj, n);
		   
  return fx;
}
*/
