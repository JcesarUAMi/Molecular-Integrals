#ifndef BOYS_H
#define BOYS_H 

#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <CL/sycl.hpp>
using namespace std;

template <typename bt> struct SplinesConstants {
  vector<bt> a;  
  vector<bt> b;  
  vector<bt> c;  
  vector<bt> d;
};

template<typename T> class Boys {
  public:
    template<typename bt> bt BoysFunction (int n, bt xj) {
      int k, ni, i, j;
      div_t m;
      bt integral, fx;

      if (xj < 1E-12) 
        xj = 0.0;
      
      if (xj == 0)
        fx = double(1.0 / (2.0*n + 1.0));
      else if (xj < 1E-4) {
        ni = n/5;
        integral = SplinesBoys<double>(ni,xj);
        ni = (ni + 1) * 5;
        fx = boysDownward<double>(integral, xj, ni, n);
      } else if (xj > 30.0) {
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

    template<typename bt> bt boysZeroDegree(bt x) {
      bt result, sqX;

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

    array<SplinesConstants<double>, 10> btv;
    array<SplinesConstants<double>, 3> btvNum;

    template<typename bt> bt boysSeriesInferior(int n, bt x) {
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

    template<typename bt> bt SplinesBoys(int li, bt xj) {
      int ind_j, val;
      bt fx, xi, sxi;
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

    template <typename bt> bt boysDownward(bt integral, bt x, int m, int p) {
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

    template <typename bt> bt boysUpward(bt x, int n) {
      T result, value;
      value = 1.0/(2.0*x);
      result = SplinesBoys<T>(9, x);
      for (int i = 50; i < n; i++) {
        result *= value * (2*i + 1);
        result -= value * exp(-x);
      }
      return result;
    }

		T sqPi = 0.5*sqrt(M_PI);
    
    template<typename bt> bt factorialDouble(int n) {
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
/*
  template<typename bt> void LectureBoysValuesTest (fstream& readFromFile, fstream& readFromFile1) {

  bt ta, tb, tc, td;
  string file = "FinalCubicSplineData.dat";
  ifstream boysFile(file.c_str(), ios::binary | ios::in);
  if (!boysFile) {
    cerr << "File couldn't be opened" << endl;
    exit (EXIT_FAILURE);
  }

  for (int i=0; i<10; i++) {
    btv[i].a.resize(300001);
    btv[i].b.resize(300001);
    btv[i].c.resize(300001);
    btv[i].d.resize(300001);
    for(int j=0; j<300001; j++) {
      boysFile.read(reinterpret_cast<char *>(&ta), sizeof(ta));
      boysFile.read(reinterpret_cast<char *>(&tb), sizeof(tb));
      boysFile.read(reinterpret_cast<char *>(&tc), sizeof(tc));
      boysFile.read(reinterpret_cast<char *>(&td), sizeof(td));
      btv[i].a[j] = ta;
      btv[i].b[j] = tb;
      btv[i].c[j] = tc;
      btv[i].d[j] = td;
    }
  }

  boysFile.close();

  string file1 = "CubicSplineDataF1num.dat";
  ifstream boysFile1(file1.c_str(), ios::binary | ios::in);
  if (!boysFile1) {
    cerr << "File couldn't be opened" << endl;
    exit (EXIT_FAILURE);
  }
  for (int i=0; i<3; i++) {
    btvNum[i].a.resize(300001);
    btvNum[i].b.resize(300001);
    btvNum[i].c.resize(300001);
    btvNum[i].d.resize(300001);
    for(int j=0; j<300001; j++) {
      boysFile1.read(reinterpret_cast<char *>(&ta), sizeof(ta));
      boysFile1.read(reinterpret_cast<char *>(&tb), sizeof(tb));
      boysFile1.read(reinterpret_cast<char *>(&tc), sizeof(tc));
      boysFile1.read(reinterpret_cast<char *>(&td), sizeof(td));
      btvNum[i].a[j] = ta;
      btvNum[i].b[j] = tb;
      btvNum[i].c[j] = tc;
      btvNum[i].d[j] = td;
    }
  }

  boysFile1.close();

}
*/

};


#endif
