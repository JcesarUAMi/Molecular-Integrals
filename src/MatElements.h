#ifndef MATRIXELE_H
#define MATRIXELE_H 

#include "Basis.h"
#include "Geometry.h"
#include "Coefficients.h"
#include "readDataFiles.h"
#include "BoysFunction.h"
#include "RysDupuisKing.h"
#include "McMurchie.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <CL/sycl.hpp>
using namespace sycl;
using namespace std;

template <typename T> 
struct atomicOrbitals {
  array<T, 3> coord;
  int ang[3];
  T expo;
  int atomCenter;
  int funcType;
  int angMom;
};

template <typename U>
struct atmInfo {
  int atmType;
  array<U, 3> coords;
};

template <typename U>
struct schwarz {
  int ind_a;
  int ind_b;
  int atmCenter1;
  int atmCenter2;
  U scr;
  U ext;
  U expo1;
  U expo2;
  U expoTot;
  U mu;
  U AB;
  U Kab;
  int angT;
  array<int, 3> angs;
  array<int, 3> ang1;
  array<int, 3> ang2;
  array<U, 3> P;
  array<U, 3> PA;
  array<U, 3> PB;
};

template <typename U> class storm {
  public:
  int factorialDouble(int n) {
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

    template <typename T> T boysZeroDegree(T x) {
      T result, sqX;

      if (x == 0.0)
        result = 1.0;
      else {
        sqX = sycl::sqrt(x);
        if (sqX >= 6.0)
          result = sqPi/sqX;
        else
          result = sycl::erf(sqX)*sqPi/sqX;
      }
      return result;
    }

  double McMurchieIntegral (int n, double redExpo, double Rpq) {

    double integral;
    double factor;
    factor = redExpo * Rpq;
    if (n == 0) {
      integral = boysZeroDegree(factor);
    } else {
      integral = mypow(-2.0*redExpo, n);
      integral *= boysSeriesInferior(n, factor);
    }

  return integral;

  }

  double mypow(double factor, int n) {

    double result;

    result = 1.0;
    if (n != 0) {
      if (n < 0) {
        n *= (-1); 
        for (int i=0; i<n; i++)
          result *= (1.0/factor);
      } else {
        for (int i=0; i<n; i++)
          result *= factor;
      }
    }

   return result; 
    

  }

  template<typename bt> bt ObaraIntegral (int n, double factor) {

  double integral;

  if (factor == 0) {
    integral = double(1.0 / (2.0*n + 1.0));
  } else if (factor >= 30.0) {
    integral = factorialDouble(2*n-1)/mypow(2.0,n+1);
    integral *= sycl::sqrt(M_PI/mypow(factor,2*n+1));
  } 

  return integral;

}

  int factorial1 (int t, int ni) {

  int result;
  result = t;

  for (int i=1; i<ni; i++)
    result *= (t-i);

  return result;

}

  int factorial (int n) {

  int result;
  result = 1;

  if (n > 1) {
    for (int i=0; i<n; i++)
      result *= (n-i);
  }

  return result;

}


  int binomialCoef (int n, int k) {

  int result;
  result = factorial(n)/(factorial(k)*factorial(n-k));

  return result;

}

int lambda (int i, int t) {

  int result;
  i *= 2;
  result = factorial(t+i) / (factorialDouble(i)*factorial(t));

  return result;
}

  int coeficient (int la, int j) {

  int result;
  if (j > la)
    result = 0;
  else
    result = factorial(la) / (factorial(j)*factorial(la-j));

  return result;
}

double aCoef (int la, int lb, int ff, int i, double Xpa, double Xpb) {

  double coef1;
  int nff;

  coef1 = 0.0;
  for (int j=0; j<=ff-2*i; j++) {
    nff = ff-2*i-j;
    coef1 += coeficient(la, j)*coeficient(lb, nff)*mypow(Xpa,j)*mypow(Xpb,nff);
  }

  return coef1;
}

double McMurchieCoefTiam (int la, int lb, int t, double p, double Xpa, double Xpb) {

  double coef1, factor, result, coef, fac;
  int ff, I;

  if (t < 0 || t > la+lb)
    result = 0.0;
  else {
    ff = la + lb - t;
    I = ff * 0.5;
    factor = mypow(2*p, -t);
    fac = 0.0;
    for (int i=0; i<=I; i++) {
      coef = lambda(i, t) * mypow(2*p, -i);
      coef1 = aCoef(la, lb, ff, i, Xpa, Xpb);
      fac += coef*coef1;
    }
    result = factor * fac;
  }

  return result;
}

  double OverlapNumericalIntegral(int la, int lb, double Xpa, double Xpb, double q) {

    int i;
    double integral = 0.0, sq;
    double x, xk[6], wk[6];
    sq = sycl::sqrt(q);
/*
  xk[0]  = -3.8897248978697819192716427472441917601481810141665;     xk[1]  = -3.0206370251208897717106793751767636460038065983103;
  xk[2]  = -2.279507080501059900187728569424341135039156593214;     xk[3]  = -1.5976826351526047967096627709045767082402869887687;
  xk[4]  = -0.94778839124016374370457813106013651367915379591455;    xk[5]  = -0.31424037625435911127661163409533712897211083832271;
  xk[6]  = 3.8897248978697819192716427472441917601481810141665;     xk[7]  = 3.0206370251208897717106793751767636460038065983103;
  xk[8]  = 2.279507080501059900187728569424341135039156593214;     xk[9]  = 1.5976826351526047967096627709045767082402869887687;
  xk[10]  = 0.94778839124016374370457813106013651367915379591455;    xk[11]  = 0.31424037625435911127661163409533712897211083832271;
  wk[0] = 2.658551684356301606023114008768679442902629821966E-7;     wk[1]  = 8.573687043587858654569063231532409401897054469197E-5;
  wk[2] = 0.00390539058462906185999438432619531097415000471400858;     wk[3]  = 0.051607985615883929991873442360613094481551034950952;
  wk[4] = 0.26049231026416112923339613976524753324486296335793;      wk[5]  = 0.57013523626247957834711348227480045173624746423056;
  wk[6] = 2.658551684356301606023114008768679442902629821966E-7;     wk[7]  = 8.573687043587858654569063231532409401897054469197E-5;
  wk[8] = 0.00390539058462906185999438432619531097415000471400858;     wk[9]  = 0.051607985615883929991873442360613094481551034950952;
  wk[10] = 0.26049231026416112923339613976524753324486296335793;      wk[11]  = 0.57013523626247957834711348227480045173624746423056;
*/

  xk[0] = -2.3506049736744922228339219870609208465378638272941;
  xk[1] = -1.335849074013696949714895282970370673972589526125;
  xk[2] = -0.43607741192761650867921594825062491076625666893183;
  xk[3] = 0.43607741192761650867921594825062491076625666893183;
  xk[4] = 1.335849074013696949714895282970370673972589526125;
  xk[5] = 2.3506049736744922228339219870609208465378638272941;
  wk[0] = 0.00453000990550884564085747256462715093251007194044;
  wk[1] = 0.15706732032285664391631156350837826891231465003571;
  wk[2] = 0.72462959522439252409191470559756717155395000608505;
  wk[3] = 0.72462959522439252409191470559756717155395000608505;
  wk[4] = 0.15706732032285664391631156350837826891231465003571;
  wk[5] = 0.00453000990550884564085747256462715093251007194044;

  for (i = 0; i < 6; i++) {
    x = xk[i]/sq;
    integral += wk[i] * mypow(x + Xpa, la) * mypow(x + Xpb, lb);
  }

  return integral * sycl::sqrt(1.0/M_PI);

}

double overlapIntegral(int la, int lb, double Xpa, double Xpb, double q) {

  double integral;
  integral = sycl::sqrt(M_PI/q);
  
  if (la == 0 && lb == 0)
    return integral;
  else
    return integral * OverlapNumericalIntegral(la, lb, Xpa, Xpb, q); 

}

double otherOverlap(int la, int lb, double Xpa, double Xpb, double q, int *Efactors) {

  double integral;
  integral = sycl::sqrt(M_PI/q);
  
  if (la == 0 && lb == 0)
    return integral;
  else
    return integral * McMurchieCoefLoop(Efactors, la, lb, 0, q, Xpa, Xpb);

}

double McMurchieCoefNew(int la, int lb, int t, double p, double Xpa, double Xpb) {

  double integral;
  double integral1, integral2;
  if (la < 0 || lb < 0)
    return 0.0;
  else if (t == 0) {
    if (la == 0 && lb == 0)
      integral = 1.0;
    else
      integral = OverlapNumericalIntegral(la, lb, Xpa, Xpb, p);
  } else if (t == 1) {
      integral = 1/(2.0*p)*(la*OverlapNumericalIntegral(la-1, lb, Xpa, Xpb, p) + lb*OverlapNumericalIntegral(la, lb-1, Xpa, Xpb, p));
  } else if (t == 2) {
    integral1 = 1/(2.0*p)*((la-1)*OverlapNumericalIntegral(la-2, lb, Xpa, Xpb, p) + lb*OverlapNumericalIntegral(la-1, lb-1, Xpa, Xpb, p));
    integral2 = 1/(2.0*p)*(la*OverlapNumericalIntegral(la-1, lb-1, Xpa, Xpb, p) + (lb-1)*OverlapNumericalIntegral(la, lb-2, Xpa, Xpb, p));
    integral = 1.0/(4.0*p)*(la * integral1 + lb * integral2);
  }

  return integral;

}

int otherFact (int la, int lb, int t) {

  int res;
  res = 1.0;
  for (int i=0; i<t; i++) {
    res *= (la-i);
  }

  if (t >= 2)
    t--;

  for (int i=0; i<=t-1; i++) {
    if (lb > i)
      res *= (lb-i);

  }

  return res;

}


int inverFact (int n, int k) {

  int result;
  result = 1.0;
  for (int i=0; i<=k; i++)
    result *= (n-i);

  return result;
}

void factors (int *Efactors, int ni, int la, int lb) {

  int nj;
  nj = ni - 1;
  Efactors[0] = inverFact(la, nj);
  Efactors[ni] = inverFact(lb, nj);
  if (ni == 2)
    Efactors[1] = la*lb;
  else if (ni == 3) {
    Efactors[1] = otherFact(la, lb, nj);
    Efactors[2] = otherFact(lb, la, nj);
  } else if (ni == 4) {
    Efactors[1] = la*lb*(la-2)*(la-1);
    Efactors[2] = la*lb*(la-1)*(lb-1);
    Efactors[3] = la*lb*(lb-1)*(lb-2);
  } else if (ni == 5) {
    Efactors[1] = la*lb*(la-1)*(la-2)*(la-3);
    Efactors[2] = la*(la-1)*lb*(la-2)*(lb-1);
    Efactors[3] = lb*(lb-1)*la*(lb-2)*(la-1);
    Efactors[4] = lb*la*(lb-1)*(lb-2)*(lb-3);
  } else if (ni == 6) {
    Efactors[1] = la*lb*(la-1)*(la-2)*(la-3)*(la-4);
    Efactors[2] = la*lb*(la-1)*(lb-1)*(la-2)*(la-3);
    Efactors[3] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2);
    Efactors[4] = lb*la*(lb-1)*(la-1)*(lb-2)*(lb-3);
    Efactors[5] = lb*la*(lb-1)*(lb-2)*(lb-3)*(lb-4);
  } else if (ni == 7) {
    Efactors[1] = la*lb*(la-1)*(la-2)*(la-3)*(la-4)*(la-5);
    Efactors[2] = la*lb*(la-1)*(lb-1)*(la-2)*(la-3)*(la-4);
    Efactors[3] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3);
    Efactors[4] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(lb-3);
    Efactors[5] = lb*la*(lb-1)*(la-1)*(lb-2)*(lb-3)*(lb-4);
    Efactors[6] = lb*la*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5);
  } else if (ni == 8) {
    Efactors[1] = la*lb*(la-1)*(la-2)*(la-3)*(la-4)*(la-5)*(la-6);
    Efactors[2] = la*lb*(la-1)*(lb-1)*(la-2)*(la-3)*(la-4)*(la-5);
    Efactors[3] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3)*(la-4);
    Efactors[4] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3)*(lb-3);
    Efactors[5] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(lb-3)*(lb-4);
    Efactors[6] = lb*la*(lb-1)*(la-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5);
    Efactors[7] = lb*la*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5)*(lb-6);
  } else if (ni == 9) {
    Efactors[1] = la*lb*(la-1)*(la-2)*(la-3)*(la-4)*(la-5)*(la-6)*(la-7);
    Efactors[2] = la*lb*(la-1)*(lb-1)*(la-2)*(la-3)*(la-4)*(la-5)*(la-6);
    Efactors[3] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3)*(la-4)*(la-5);
    Efactors[4] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3)*(lb-3)*(la-4);
    Efactors[5] = lb*la*(la-1)*(lb-1)*(la-2)*(la-3)*(lb-2)*(lb-3)*(lb-4);
    Efactors[6] = lb*la*(lb-1)*(la-1)*(la-2)*(lb-2)*(lb-3)*(lb-4)*(lb-5);
    Efactors[7] = lb*la*(la-1)*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5)*(lb-6);
    Efactors[8] = lb*la*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5)*(lb-6)*(lb-7);
  } else if (ni == 10) {
    Efactors[1] = la*lb*(la-1)*(la-2)*(la-3)*(la-4)*(la-5)*(la-6)*(la-7)*(la-8);
    Efactors[2] = la*lb*(la-1)*(lb-1)*(la-2)*(la-3)*(la-4)*(la-5)*(la-6)*(la-7);
    Efactors[3] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3)*(la-4)*(la-5)*(la-6);
    Efactors[4] = la*lb*(la-1)*(lb-1)*(la-2)*(lb-2)*(la-3)*(lb-3)*(la-4)*(la-5);
    Efactors[5] = lb*la*(la-1)*(lb-1)*(la-2)*(la-3)*(la-4)*(lb-2)*(lb-3)*(lb-4);
    Efactors[6] = lb*la*(lb-1)*(la-1)*(la-2)*(la-3)*(lb-2)*(lb-3)*(lb-4)*(lb-5);
    Efactors[7] = lb*la*(la-1)*(la-2)*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5)*(lb-6);
    Efactors[8] = lb*la*(la-1)*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5)*(lb-6)*(lb-7);
    Efactors[9] = lb*la*(lb-1)*(lb-2)*(lb-3)*(lb-4)*(lb-5)*(lb-6)*(lb-7)*(lb-8);

  }

}

  double McMurchieCoefLoop (int *Efactors, int la, int lb, int t, double p, double Xpa, double Xpb) {

  int nT, ni, nj, value, i, j;
  double integral;
  double denominator;

  if (t < 0 || t > la+lb)
    return 0.0;
  else if (la == 0 && lb > 0) {
      integral = mypow(2.0*p,-t)*binomialCoef(lb, t) * OverlapNumericalIntegral(0, lb-t, Xpa, Xpb, p);
  } else if (lb == 0 && la > 0) {
      integral = mypow(2.0*p,-t)*binomialCoef(la, t) * OverlapNumericalIntegral(la-t, 0, Xpa, Xpb, p);
  } else if (t <= 2) {
    integral = McMurchieCoefNew(la, lb, t, p, Xpa, Xpb);
  } else {
    ni = t - 2;
    denominator = 1.0/(mypow(2.0*p, ni) * factorial1(t, ni));
    factors(Efactors, ni, la, lb);
    integral = 0.0;
    for (int k=0; k<=ni; k++) {
      value = binomialCoef(ni, k);
      j = lb - k;
      i = la - (ni - k);
      integral += value * Efactors[k] * McMurchieCoefNew(i, j, 2, p, Xpa, Xpb);
    }
    integral *= denominator;
  }

  return integral;
}

void DerivativeBoysOneCenter(double *boysD, int *boysN, double p, int n, int& ni) {

float dif;
ni = 1;

if (n == 0) {
  boysD[0] = 1.0;
  boysN[0] = 0;
} else {
  dif = float(n/2.0) - n/2;
  if (dif != 0) {
  boysD[0] = 0.0;
  boysN[0] = 0;
} else {
  if (n == 2) {
    boysD[0] = -2.0*p;
    boysN[0] = 1;
  } else if (n == 4) {
    boysD[0] = 12.0*p*p;
    boysN[0] = 2;
  } else if (n == 6) {
    boysD[0] = -120.0*p*p*p;
    boysN[0] = 3;
  } else if (n == 8) {
    boysD[0] = 1680.0*p*p*p*p;
    boysN[0] = 4;
  } else if (n == 10) {
    boysD[0] = -30240.0*mypow(p,5);
    boysN[0] = 5;
  } else if (n == 12) {
    boysD[0] = 665280.0*mypow(p,6);
    boysN[0] = 6;
  } else if (n == 14) {
    boysD[0] = -17297280 * mypow(p,7);
    boysN[0] = 7;
  } else if (n == 16) {
    boysD[0] = 518918400 * mypow(p,8);
    boysN[0] = 8;
  }
  
  }
}

}

 void DerivativeBoys (double *boysD, int *boysN, double p, double Xpq, int n, int& ni) {

  double preFac;
  preFac = -2.0*p*Xpq;

  if (n == 0) {
    ni = 1;
    boysD[0] = 1.0;
    boysN[0] = 0;
  } else if (n == 1) {
    ni = 1;
    boysD[0] = preFac;
    boysN[0] = 1;
  } else if (n == 2) {
    ni = 2;
    boysD[0] = -2.0*p;
    boysN[0] = 1;
    boysD[1] = preFac*preFac;
    boysN[1] = 2;
  } else if (n == 3) {
    ni = 2;
    boysD[0] = 12.0*p*p*Xpq;
    boysN[0] = 2;
    boysD[1] = preFac*preFac*preFac;
    boysN[1] = 3;
  } else if (n == 4) {
    ni = 3;
    boysD[0] = 12.0*p*p;
    boysN[0] = 2;
    boysD[1] = -48.0*p*p*p*Xpq*Xpq;
    boysN[1] = 3;
    boysD[2] = preFac*preFac*preFac*preFac;
    boysN[2] = 4;
  } else if (n == 5) {
    ni = 3;
    boysD[0] = -120.0*p*p*p*Xpq;
    boysN[0] = 3;
    boysD[1] = 160.0*p*p*p*p*Xpq*Xpq*Xpq;
    boysN[1] = 4;
    boysD[2] = preFac*preFac*preFac*preFac*preFac;
    boysN[2] = 5;
  } else if (n == 6) {
    ni = 4;
    boysD[0] = -120.0*p*p*p;
    boysD[1] = 720.0*p*p*p*p*Xpq*Xpq;
    boysD[2] = -480.0*p*p*p*p*p*Xpq*Xpq*Xpq*Xpq;
    boysD[3] = mypow(preFac,6);
    boysN[0] = 3;
    boysN[1] = 4;
    boysN[2] = 5;
    boysN[3] = 6;
   } else if (n == 7) {
    ni = 4;
    boysD[0] = 1680.0*p*p*p*p*Xpq;
    boysD[1] = -3360.0*p*p*p*p*p*Xpq*Xpq*Xpq;
    boysD[2] = 1344.0*mypow(p,6)*mypow(Xpq,5);
    boysD[3] = mypow(preFac,7);
    boysN[0] = 4;
    boysN[1] = 5;
    boysN[2] = 6;
    boysN[3] = 7;
   } else if (n == 8) {
    ni = 5;
    boysD[0] = 1680.0*p*p*p*p;
    boysD[1] = -13440*mypow(p,5)*Xpq*Xpq;
    boysD[2] = 13440.0*mypow(p,6)*Xpq*Xpq*Xpq*Xpq;
    boysD[3] = -3584.0*mypow(p,7)*mypow(Xpq,6);
    boysD[4] = mypow(preFac,8);
    boysN[0] = 4;
    boysN[1] = 5;
    boysN[2] = 6;
    boysN[3] = 7;
    boysN[4] = 8;
   } else if (n == 9) {
    ni = 5;
    boysD[0] = -30240.0*mypow(p,5)*Xpq;
    boysD[1] = 80640.0*mypow(p,6)*Xpq*Xpq*Xpq;
    boysD[2] = -48384.0*mypow(p,7)*mypow(Xpq,5);
    boysD[3] = 9216.0*mypow(p,8)*mypow(Xpq,7);
    boysD[4] = mypow(preFac,9);
    boysN[0] = 5;
    boysN[1] = 6;
    boysN[2] = 7;
    boysN[3] = 8;
    boysN[4] = 9;
   } else if (n == 10) {
   ni = 6;
    boysD[0] = -30240.0*mypow(p,5);
    boysD[1] = 302400.0*mypow(p,6)*Xpq*Xpq;
    boysD[2] = -403200.0*mypow(p,7)*Xpq*Xpq*Xpq*Xpq;
    boysD[3] = 161280.0*mypow(p,8)*mypow(Xpq,6);
    boysD[4] = -23040.0*mypow(p,9)*mypow(Xpq,8);
    boysD[5] = mypow(preFac,10);
    boysN[0] = 5;
    boysN[1] = 6;
    boysN[2] = 7;
    boysN[3] = 8;
    boysN[4] = 9;
    boysN[5] = 10;
  } else if (n == 11) {
    ni = 6;
    boysD[0] = 665280.0*mypow(p,6)*Xpq;
    boysD[1] = -2217600.0*mypow(p,7)*Xpq*Xpq*Xpq;
    boysD[2] = 1774080.0*mypow(p,8)*mypow(Xpq,5);
    boysD[3] = -506880.0*mypow(p,9)*mypow(Xpq,7);
    boysD[4] = 56320.0*mypow(p,10)*mypow(Xpq,9);
    boysD[5] = mypow(preFac, 11);
    boysN[0] = 6;
    boysN[1] = 7;
    boysN[2] = 8;
    boysN[3] = 9;
    boysN[4] = 10;
    boysN[5] = 11;
  } else if (n == 12) {
    ni = 7;
    boysD[0] = 665280.0*mypow(p,6);
    boysD[1] = -7983360.0*mypow(p,7)*Xpq*Xpq;
    boysD[2] = 13305600.0*mypow(p,8)*Xpq*Xpq*Xpq*Xpq;
    boysD[3] = -7096320.0*mypow(p,9)*mypow(Xpq,6);
    boysD[4] = 1520640.0*mypow(p,10)*mypow(Xpq,8);
    boysD[5] = -135168.0*mypow(p,11)*mypow(Xpq,10);
    boysD[6] = mypow(preFac,12);
    boysN[0] = 6;
    boysN[1] = 7;
    boysN[2] = 8;
    boysN[3] = 9;
    boysN[4] = 10;
    boysN[5] = 11;
    boysN[6] = 12;
  } else if (n == 13) {
    ni = 7;
    boysD[0] = -17297280.0*mypow(p,7)*Xpq;
    boysD[1] = 69189120.0*mypow(p,8)*Xpq*Xpq*Xpq;
    boysD[2] = -69189120.0*mypow(p,9)*mypow(Xpq,5);
    boysD[3] = 26357760.0*mypow(p,10)*mypow(Xpq,7);
    boysD[4] = -4392960.0*mypow(p,11)*mypow(Xpq,9);
    boysD[5] = 319488.0*mypow(p,12)*mypow(Xpq,11);
    boysD[6] = mypow(preFac,13);
    boysN[0] = 7;
    boysN[1] = 8;
    boysN[2] = 9;
    boysN[3] = 10;
    boysN[4] = 11;
    boysN[5] = 12;
    boysN[6] = 13;
  } else if (n == 14) {
    ni = 8;
    boysD[0] = -17297280 * mypow(p,7);
    boysD[1] = 242161920 * mypow(p,8)*Xpq*Xpq;
    boysD[2] = -484323840 * mypow(p,9)*mypow(Xpq,4);
    boysD[3] = 322882560 * mypow(p,10)*mypow(Xpq,6);
    boysD[4] = -92252160 * mypow(p,11)*mypow(Xpq,8);
    boysD[5] = 12300288 * mypow(p,12)*mypow(Xpq,10);
    boysD[6] = -745472 * mypow(p,13)*mypow(Xpq,12);
    boysD[7] = mypow(preFac,14);
    boysN[0] = 7;
    boysN[1] = 8;
    boysN[2] = 9;
    boysN[3] = 10;
    boysN[4] = 11;
    boysN[5] = 12;
    boysN[6] = 13;
    boysN[7] = 14;
  } else if (n == 15) {
    ni = 8;
    boysD[0] = 518918400 * mypow(p,8)*Xpq;
    boysD[1] = -2421619200 * mypow(p,9)*Xpq*Xpq*Xpq;
    boysD[2] = 2905943040 * mypow(p,10)*mypow(Xpq,5);
    boysD[3] = -1383782400 * mypow(p,11)*mypow(Xpq,7);
    boysD[4] = 307507200 * mypow(p,12)*mypow(Xpq,9);
    boysD[5] = -33546240 * mypow(p,13)*mypow(Xpq,11);
    boysD[6] = 1720320 * mypow(p,14)*mypow(Xpq,13);
    boysD[7] = mypow(preFac,15);
    boysN[0] = 8;
    boysN[1] = 9;
    boysN[2] = 10;
    boysN[3] = 11;
    boysN[4] = 12;
    boysN[5] = 13;
    boysN[6] = 14;
    boysN[7] = 15;
  } else if ( n == 16) {
    ni = 9;
    boysD[0] = 518918400 * mypow(p,8);
    boysD[1] = -8302694400 * mypow(p,9)*Xpq*Xpq;
    boysD[2] = 19372953600 * mypow(p,10)*mypow(Xpq,4);
    boysD[3] = -15498362880 * mypow(p,11)*mypow(Xpq,6);
    boysD[4] = 5535129600 * mypow(p,12)*mypow(Xpq,8);
    boysD[5] = -984023040 * mypow(p,13)*mypow(Xpq,10);
    boysD[6] = 89456640 * mypow(p,14)*mypow(Xpq,12);
    boysD[7] = -3932160 * mypow(p,15)*mypow(Xpq,14);
    boysD[8] = mypow(preFac,16);
    boysN[0] = 8;
    boysN[1] = 9;
    boysN[2] = 10;
    boysN[3] = 11;
    boysN[4] = 12;
    boysN[5] = 13;
    boysN[6] = 14;
    boysN[7] = 15;
    boysN[8] = 16;
  } else if ( n == 17) {
    ni = 9;
    boysD[0] =  -17643225600 * mypow(p,9) * Xpq;
    boysD[1] = 94097203200 * mypow(p,10)*Xpq*Xpq*Xpq;
    boysD[2] = -131736084480 * mypow(p,11)*mypow(Xpq,5);
    boysD[3] = 75277762560 * mypow(p,12)*mypow(Xpq,7);
    boysD[4] =  -20910489600 * mypow(p,13)*mypow(Xpq,9);
    boysD[5] =  3041525760 * mypow(p,14)*mypow(Xpq,11);
    boysD[6] = -233963520 * mypow(p,15)*mypow(Xpq,13);
    boysD[7] = 8912896 * mypow(p,16)*mypow(Xpq,15);
    boysD[8] = mypow(preFac,17);
    boysN[0] = 9;
    boysN[1] = 10;
    boysN[2] = 11;
    boysN[3] = 12;
    boysN[4] = 13;
    boysN[5] = 14;
    boysN[6] = 15;
    boysN[7] = 16;
    boysN[8] = 17;
  } else if ( n == 18) {
    ni = 10;
    boysD[0] =  -17643225600 * mypow(p,9);
    boysD[1] = 317578060800 * mypow(p,10)*Xpq*Xpq;
    boysD[2] = - 846874828800* mypow(p,11)*mypow(Xpq,4);
    boysD[3] = 790416506880 * mypow(p,12)*mypow(Xpq,6);
    boysD[4] = -338749931520 * mypow(p,13)*mypow(Xpq,8);
    boysD[5] = 75277762560 * mypow(p,14)*mypow(Xpq,10);
    boysD[6] = -9124577280* mypow(p,15)*mypow(Xpq,12);
    boysD[7] = 601620480 * mypow(p,16)*mypow(Xpq,14);
    boysD[8] = -20054016 * mypow(p,17)*mypow(Xpq,16);
    boysD[9] = mypow(preFac,18);
    boysN[0] = 9;
    boysN[1] = 10;
    boysN[2] = 11;
    boysN[3] = 12;
    boysN[4] = 13;
    boysN[5] = 14;
    boysN[6] = 15;
    boysN[7] = 16;
    boysN[8] = 17;
    boysN[9] = 18;
  } else if ( n == 19) {
    ni = 10;
    boysD[0] = 670442572800 * mypow(p,10) * Xpq;
    boysD[1] = -4022655436800 * mypow(p,11)*Xpq*Xpq*Xpq;
    boysD[2] = 6436248698880 * mypow(p,12)*mypow(Xpq,5);
    boysD[3] = -4290832465920* mypow(p,13)*mypow(Xpq,7);
    boysD[4] = 1430277488640 * mypow(p,14)*mypow(Xpq,9);
    boysD[5] = -260050452480* mypow(p,15)*mypow(Xpq,11);
    boysD[6] = 26671841280* mypow(p,16)*mypow(Xpq,13);
    boysD[7] = -1524105216* mypow(p,17)*mypow(Xpq,15);
    boysD[8] = 44826624* mypow(p,18)*mypow(Xpq,17);
    boysD[9] = mypow(preFac,19);
    boysN[0] = 10;
    boysN[1] = 11;
    boysN[2] = 12;
    boysN[3] = 13;
    boysN[4] = 14;
    boysN[5] = 15;
    boysN[6] = 16;
    boysN[7] = 17;
    boysN[8] = 18;
    boysN[9] = 19;
  } else if ( n == 20) {
    ni = 11;
    boysD[0] = 670442572800*mypow(p,10);
    boysD[1] = -13408851456000* mypow(p,11)*Xpq*Xpq;
    boysD[2] = 40226554368000* mypow(p,12)*mypow(Xpq,4);
    boysD[3] = -42908324659200* mypow(p,13)*mypow(Xpq,6);
    boysD[4] = 21454162329600* mypow(p,14)*mypow(Xpq,8);
    boysD[5] = -5721109954560* mypow(p,15)*mypow(Xpq,10);
    boysD[6] = 866834841600* mypow(p,16)*mypow(Xpq,12);
    boysD[7] = -76205260800* mypow(p,17)*mypow(Xpq,14);
    boysD[8] = 3810263040* mypow(p,18)*mypow(Xpq,16);
    boysD[9] = -99614720* mypow(p,19)*mypow(Xpq,18);
    boysD[10] = mypow(preFac,20);
    boysN[0] = 10;
    boysN[1] = 11;
    boysN[2] = 12;
    boysN[3] = 13;
    boysN[4] = 14;
    boysN[5] = 15;
    boysN[6] = 16;
    boysN[7] = 17;
    boysN[8] = 18;
    boysN[9] = 19;
    boysN[10] = 20;
  }

 
}

  double matrixRes (double redExpo, double Xpq, double Ypq, double Zpq, double *boysDx, double *boysDy, double *boysDz, int *boysNx, int *boysNy, int *boysNz) {

  double other, res;
  int n;
  int nx, ny, nz;

  res = 0.0;
  for (int mx=0; mx<nx; mx++) 
    for (int my=0; my<ny; my++) 
      for (int mz=0; mz<nz; mz++) {
        other = boysDx[mx]*boysDy[my]*boysDz[mz];
        n = boysNx[mx]+boysNy[my]+boysNz[mz];
        res += other;
      }

  return res;
}


  template<typename bt> bt boysSeriesInferior(int n, bt x) {
      U numerator, denominator, result;
      numerator = tgamma(n+0.5);
      result = 0.0;
      for (int i=0; i<50; i++) {
        denominator = tgamma(n+i+1.5);
        result += mypow(x, i)*numerator/denominator;
      }
      result *= 0.5 * sycl::exp(-x);

      return result;
    }

    template<typename bt> bt SplinesBoys(bt xi, bt vala, bt valb, bt valc, bt vald) {

      bt fx, sxi, fac1;

      fac1 = valc + vald*xi;

      fx = vala;
      fx += xi*(valb + xi*fac1);

      return fx;
    }

  inline double boysDownward(double x, double boys, int m, int p, double *boysT) {
    U result;
    int n, j;
    j = m - p;
    result = boys;
    boysT[m] = result;
    for (int i = 0; i < j; i++) {
      n = m - (i+1);
      result *= 2 * x;
      result += sycl::exp(-x);
      result /= (2 * n + 1);
      boysT[n] = result;
    }

    return result;
  
  }

    template <typename bt> bt boysUpward(bt x, int n, bt integral) {
      U result, value;
      value = 1.0/(2.0*x);
      result = integral;
      for (int i = 50; i < n; i++) {
        result *= value * (2*i + 1);
        result -= value * sycl::exp(-x);
      }
      return result;
    }

    U sqPi = 0.5*sycl::sqrt(M_PI);

		U pi52 = 2.0*sycl::sqrt(M_PI)*M_PI*M_PI;

    U tolerance = 1E-8;

    template<typename T> T OneCenterNewIntegral (double *boysT, double *boysDx, double *boysDy, double *boysDz, int *boysNx, int *boysNy, int *boysNz, int *EfactorsX, int *EfactorsY, int *EfactorsZ, int *ang12, int *ang34, int *ang1, int *ang2, int *ang3, int *ang4, int totAngMom, double etaCou, double expo12, double expo34) {

  double bielecTot, coef1Ex, coef1Ey, coef1Ez, coef2Ex, coef2Ey, coef2Ez, res, integralBie;
  int nx, ny, nz, ni;

  boysT[0] = 1.0;
  for (int i=1; i<=totAngMom; i++)
    boysT[i] = ObaraIntegral<double>(i, 0.0);

  bielecTot = 0.0;
  for (int t=0; t<=ang12[0]; t++)
    for (int u=0; u<=ang12[1]; u++)
      for (int v=0; v<=ang12[2]; v++) {
        coef1Ex = McMurchieCoefTiam(ang1[0], ang2[0], t, expo12, 0.0, 0.0);
        coef1Ey = McMurchieCoefTiam(ang1[1], ang2[1], u, expo12, 0.0, 0.0);
        coef1Ez = McMurchieCoefTiam(ang1[2], ang2[2], v, expo12, 0.0, 0.0);
        for (int ta=0; ta<=ang34[0]; ta++)
          for (int ua=0; ua<=ang34[1]; ua++)
            for (int va=0; va<=ang34[2]; va++) {
              coef2Ex = McMurchieCoefTiam(ang3[0], ang4[0], ta, expo34, 0.0, 0.0);
              coef2Ey = McMurchieCoefTiam(ang3[1], ang4[1], ua, expo34, 0.0, 0.0);
              coef2Ez = McMurchieCoefTiam(ang3[2], ang4[2], va, expo34, 0.0, 0.0);
              DerivativeBoysOneCenter(boysDx, boysNx, etaCou, t+ta, nx);
              DerivativeBoysOneCenter(boysDy, boysNy, etaCou, u+ua, ny);
              DerivativeBoysOneCenter(boysDz, boysNz, etaCou, v+va, nz);

              integralBie = 0.0;
              for (int mx=0; mx<nx; mx++)
                for (int my=0; my<ny; my++)
                  for (int mz=0; mz<nz; mz++) {
                    res = boysDx[mx]*boysDy[my]*boysDz[mz];
                    ni = boysNx[mx]+boysNy[my]+boysNz[mz];
                    integralBie += boysT[ni]*res;
                  }

              
              bielecTot += mypow(-1,ta+va+ua)*coef1Ex * coef1Ey * coef1Ez * coef2Ex * coef2Ey * coef2Ez * integralBie;

          }

    }

  return bielecTot;

}

    double NewIntegral (double *boysT, double *boysDx, double *boysDy, double *boysDz, int *boysNx, int *boysNy, int *boysNz, int *EfactorsX, int *EfactorsY, int *EfactorsZ, int *ang12, int *ang34, int *ang1, int *ang2, int *ang3, int *ang4, double *PQ, double *PA, double *PB, double *QC, double *QD, double etaCou, double expo12, double expo34) {

  double bielecTot, coef1Ex, coef1Ey, coef1Ez, coef2Ex, coef2Ey, coef2Ez, res, integralBie;
  int nx, ny, nz, ni;

  bielecTot = 0.0;
  for (int t=0; t<=ang12[0]; t++)
    for (int u=0; u<=ang12[1]; u++)
      for (int v=0; v<=ang12[2]; v++) {
        coef1Ex = McMurchieCoefTiam(ang1[0], ang2[0], t, expo12, PA[0], PB[0]);
        coef1Ey = McMurchieCoefTiam(ang1[1], ang2[1], u, expo12, PA[1], PB[1]);
        coef1Ez = McMurchieCoefTiam(ang1[2], ang2[2], v, expo12, PA[2], PB[2]);
        for (int ta=0; ta<=ang34[0]; ta++)
          for (int ua=0; ua<=ang34[1]; ua++)
            for (int va=0; va<=ang34[2]; va++) {
              coef2Ex = McMurchieCoefTiam(ang3[0], ang4[0], ta, expo34, QC[0], QD[0]);
              coef2Ey = McMurchieCoefTiam(ang3[1], ang4[1], ua, expo34, QC[1], QD[1]);
              coef2Ez = McMurchieCoefTiam(ang3[2], ang4[2], va, expo34, QC[2], QD[2]);
              DerivativeBoys(boysDx, boysNx, etaCou, PQ[0], t+ta, nx);
              DerivativeBoys(boysDy, boysNy, etaCou, PQ[1], u+ua, ny);
              DerivativeBoys(boysDz, boysNz, etaCou, PQ[2], v+va, nz);

              integralBie = 0.0;
              for (int mx=0; mx<nx; mx++)
                for (int my=0; my<ny; my++)
                  for (int mz=0; mz<nz; mz++) {
                    res = boysDx[mx]*boysDy[my]*boysDz[mz];
                    ni = boysNx[mx]+boysNy[my]+boysNz[mz];
                    integralBie += boysT[ni]*res;
                  }
              bielecTot += mypow(-1,ta+va+ua)  * coef1Ex *  coef1Ey * coef1Ez  * coef2Ex * coef2Ey * coef2Ez * integralBie;

          }

    }

  return bielecTot;

}

  int checkGroupSize (int n) {
  
    int i;

    i = 256;
    if (i >= n)
      i = 128;
    else if (i >= n)
      i = 64;
    else if (i >= n)
      i = 32;
    else if (i >= n)
      i = 16;
    else if (i >= n)
      i = 8;

    while (n%i != 0) {
      i--;
    }

    return i;
  }

};

template <typename bt>
struct oneElectro {
  
  bt kE;
  bt oE;
  bt nE;
  int ind_a;
  int ind_b;
};


template <typename bt> 
struct bielectro {

  double val;
  int ind_a;
  int ind_b;
  int ind_c;
  int ind_d;
};

class MatElements {
  private:
    readDataInput readData;
    McMD MD;
		Basis basis;
    Boys<double> bF;
    storm<double> ghF;
    vector<atmInfo<double>> atoms;
		DKR dkr;
  public:
    vector<double> btvNuma;
    vector<double> btvNumb;
    vector<double> btvNumc;
    vector<double> btvNumd;
    vector<double>btva;
    vector<double>btvb;
    vector<double>btvc;
    vector<double>btvd;
    void fourIndexTransf();
    void OrbitalData(int);
    vector<atomicOrbitals<double>> aO;
    vector<double> oE;
    vector<double> kE;
    vector<double> nE;
    vector<double> cE;
    vector<bielectro<double>> bETS;
    vector<bielectro<double>> bE;
    vector<bielectro<double>> bE1;
    vector<bielectro<double>> bE2;
    vector<bielectro<double>> bE3;
    vector<bielectro<double>> bE4;
    vector<bielectro<double>> bE5;
    vector<bielectro<double>> bE6;
    vector<bielectro<double>> bE7;
    vector<oneElectro<double>> oEi;
    vector<oneElectro<double>> oETS;
    vector<schwarz<double>> schF;
    vector<schwarz<double>> sch;
    void passingElements(int, vector<bielectro<double>>, vector<double>&);
    void Data(fstream&, fstream&, fstream&, fstream&, fstream&);
    void matrixElements(fstream&, fstream&);
    void matrixBuilding();
    void overlapMatrix(int);
    vector<double> oneElecCoulombEnergy;
    int occOrb;
    vector<double> kineticEnergy;
    vector<double> nucleousElecEnergy;
    vector<double> repulsionElecEnergy;
		double eN;
    double totalKineticEnergy;
    double totalNucleousElecEnergy;;
		double pi52 = 2.0*sycl::sqrt(M_PI)*M_PI*M_PI;
    double totalCoulombicEnergy;
    double totalExchangeEnergy;
    vector<double> coulombicEnergy;
    vector<double> exchangeEnergy;
    template <auto query, typename T>
    void do_query( const T& obj_to_query, const std::string& name, int indent=4) {
      cout << std::string(indent, ' ') << name << " is '" << obj_to_query.template get_info<query>() << "\n";
    }
    void detectDevice ();
    const double icef = 4.05223724387138872771174646914005279541015625; //Inverse Complementary error function of 1E-8 

    void IntegralEvaluation (int, fstream&, fstream&); 
    
  template<typename bt> void LectureBoysValuesTest (fstream& readFromFile, fstream& readFromFile1) {

  bt ta, tb, tc, td;
  string file = "FinalCubicSplineData.dat";
  ifstream boysFile(file.c_str(), ios::binary | ios::in);
  if (!boysFile) {
    cerr << "File couldn't be opened" << std::endl;
    exit (EXIT_FAILURE);
  }

  btva.resize(3000010);
  btvb.resize(3000010);
  btvc.resize(3000010);
  btvd.resize(3000010);
  btvNuma.resize(900003);
  btvNumb.resize(900003);
  btvNumc.resize(900003);
  btvNumd.resize(900003);
  
  int u=0;
  for (int i=0; i<10; i++) {
    for(int j=0; j<300001; j++) {
      boysFile.read(reinterpret_cast<char *>(&ta), sizeof(ta));
      boysFile.read(reinterpret_cast<char *>(&tb), sizeof(tb));
      boysFile.read(reinterpret_cast<char *>(&tc), sizeof(tc));
      boysFile.read(reinterpret_cast<char *>(&td), sizeof(td));
      btva[u] = ta;
      btvb[u] = tb;
      btvc[u] = tc;
      btvd[u] = td;
      u++;
    }
  }

  boysFile.close();

  string file1 = "CubicSplineDataF1num.dat";
  ifstream boysFile1(file1.c_str(), ios::binary | ios::in);
  if (!boysFile1) {
    cerr << "File couldn't be opened" << std::endl;
    exit (EXIT_FAILURE);
  } 
  u = 0;
  for (int i=0; i<3; i++) {
    for(int j=0; j<300001; j++) {
      boysFile1.read(reinterpret_cast<char *>(&ta), sizeof(ta));
      boysFile1.read(reinterpret_cast<char *>(&tb), sizeof(tb));
      boysFile1.read(reinterpret_cast<char *>(&tc), sizeof(tc));
      boysFile1.read(reinterpret_cast<char *>(&td), sizeof(td));
      btvNuma[u] = ta;
      btvNumb[u] = tb;
      btvNumc[u] = tc;
      btvNumd[u] = td;
  //    cout << i << "   " << j << "   " << btvNuma[i][j] << "   " << btvNumb[i][j] << "   " << btvNumc[i][j] << "   " << btvNumd[i][j] << std::endl;
      u++;
    }
  }

  boysFile1.close();

}


};

#endif
