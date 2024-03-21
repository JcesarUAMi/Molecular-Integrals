#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <CL/sycl.hpp>
#include "RysDupuisKing.h"
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
using namespace std;


float DKR::funcQ (int i, int n, float r)  {

  if (n == 0)
    return 1.0;
  else if (n == 1)
    return r;
   else
    return r * funcQ(i, n-1, r) - var.B[n-1] * funcQ(i, n-2, r);

}

float DKR::funcU (float P, float J, int n) {

  if (n == 0)
    return sqrt(M_PI / l1) * erf(J);
  else
    return (n - 1)/(2*l1) * funcU(P, J, n-2) - P;
}

void DKR::funcG (int n) {

  float y, r;
  var.G[n] = 0.0;
  for (int i = 0; i < 50; i++) {
    r = leg.r[i];
    y = funcQ(i,n,r) * funcQ(i,n,r);
    var.G[n] += leg.E[i] * y;
  }
}


void DKR::funcB (int p) {

  int n;
  for (n = p; n < Ns; n++) {
    funcG(n);
    var.B[n] = var.G[n] / var.G[n-1];
  }
}

void DKR::funcJ () {
  int h;
  for (int n = 0; n < Ns; n++) {
    var.J[n] = 0.;
    h = n + 1;
    var.J[n] = sqrt(var.B[h]);
    if (n == Ns - 1)
      var.J[n] = 0.;
  }
}

float DKR::pythag (float a, float b) {

  float result;
  result = sqrt(a*a + b*b);

  return result;
}

void DKR::funcP () {
  float x, r;
  for (int n = 0; n < Ns; n++) {
    var.wi[n] = 0.0;
    x = 0.;
    for (int i = 0; i < Ns; i++) {
      r = var.R[n];
      x += funcQ(n, i, r) * funcQ(n, i, r) / var.G[i];
    }
    var.wi[n] = 1. / x;
    var.R[n] *= var.R[n];
    }
}

void DKR::shellDstev () {
  int iter, m, i, l;
  float s, p, dd, c, b;
  float g, r, f;
  for(l = 0; l < Ns; l++) {
    iter= 0;
    do {
      for (m = l; m < Ns-1; m++) {
        dd = fabs(var.R[m])+fabs(var.R[m+1]);
        if ((float)(fabs(var.J[m])+dd) == dd) break;
      }
      if(m != l) {
        if (iter++ == 30)
          cout << "try again\n";
        g = (var.R[l+1]-var.R[l])/(2.0*var.J[l]);
        r = pythag(g,1.0);
        g = var.R[m]-var.R[l]+var.J[l]/(g+SIGN(r,g));
        s = c = 1.0;
        p = 0.0;
        for (i = m-1; i >= l; i--) {
          f = s * var.J[i];
          b = c * var.J[i];
          var.J[i+1] = (r = pythag(f,g));
          if (r == 0.0) {
            var.R[i+1] -= p;
            var.J[m] = 0.0;
            break;
          }
          s = f/r;
          c = g/r;
          g = var.R[i+1] - p;
          r =(var.R[i] - g) * s + 2.0 * c * b;
          var.R[i+1] = g + (p = s*r);
          g = c * r - b;
        }
        if (r == 0.0 && i >= l)
          continue;
        var.R[l] -= p;
        var.J[l] = g;
        var.J[m] = 0.0;
      }
    } while (m != l);
  }
}

double DKR::DKRrec (int t, double Xpq, double eta, float s) {

  double integral;

  if (t == 0)
    integral =  1.0;
  else if (t == 1) {
    integral = -2.0*eta*s*Xpq;
  } else
    integral = -2.0*eta*s*(Xpq*DKRrec(t-1, Xpq, eta, s) + (t-1)*DKRrec(t-2, Xpq, eta, s));

  return integral;
}


void DKR::evalWaR() {

  float y, B123, J, P;
  int h;
	var.B.clear();
	var.G.clear();
	var.J.clear();
	var.R.clear();
	var.wi.clear();

  var.B.resize(Ns);
  var.G.resize(Ns);
  var.J.resize(Ns);
  var.R.resize(Ns);
  var.wi.resize(Ns);

		
 	if (l1 < 1.0) {
  	for (int i=0; i<50; i++) {
      y = leg.r[i] * leg.r[i];
      leg.E[i] = 2.0 * leg.w[i] * exp(-l1 * y);
      var.G[0] += leg.E[i];
      var.B[0] += leg.E[i];
    }
    funcB(1);
  } else {
    J = sqrt(l1);
    P = (1.0/l1) * exp(-l1);
    var.B[0] = sqrt(M_PI / l1) * erf(J);
    var.G[0] = funcU(P, J, 0);

    var.G[1] = funcU(P, J, 2);
    var.B[1] = var.G[1] / var.G[0];

    var.G[2] = funcU(P, J, 4) - 2.0 * var.B[1] * funcU(P, J, 2) + var.B[1] * var.B[1] * funcU(P, J, 0);
    var.B[2] = var.G[2] / var.G[1];

    var.G[3] = funcU(P, J, 6) - 2.0 * (var.B[1] + var.B[2]) * funcU(P, J, 4) + (var.B[2] + var.B[1]) * (var.B[2] + var.B[1]) * funcU(P, J, 2);
    var.B[3] = var.G[3] / var.G[2];

//    B123 = var.B[1] + var.B[2] + var.B[3];
//    var.G[4] = funcU(B0, y, 8) - 2 * B123 * funcU(B0, y, 6) + (2 * var.B[1] * var.B[3] + B123 * B123) * funcU(B0, y, 4) - 2 * B123 * var.B[1] * var.B[3] * funcU(B0, y, 2) + var.B[1] * var.B[1] * var.B[3] * var.B[3] * funcU(B0, y, 0);
//    var.B[4] = var.G[4] / var.G[3];

    if (Ns > 3) {
      for (int i=0; i<50; i++) {
        y = leg.r[i] * leg.r[i];
        leg.E[i] = 2.0 * leg.w[i] * exp(-l1 * y);
      }
      funcB(4);
    }
  }

  for (int n1 = 0; n1 < Ns; n1++) {
    var.J[n1] = 0.0;
    var.R[n1] = 0.0;
    h = n1 + 1;
    var.J[n1] = sqrt(var.B[h]);
    if (n1 == Ns - 1)
      var.J[n1] = 0.0;
  }

  shellDstev();

  funcP();

}
/*
double DKR::recurDKR (int i, int j, int k, int l, double Xpa, double Xpb, double Xqc, double Xqd, double Xpq, double eta, double p, double q, double s) {

  double integral;

  if (i < 0 || j < 0 || k < 0 || l < 0)
    integral = 0;
  else {
    if (i == 0 && j == 0 && k == 0 && l == 0)
      integral = 1.0;
    else if (j == 0 && k == 0 && l == 0) {
      if (i == 1)  
        integral = (Xpa - (eta/p)*Xpq*s)*recurDKR(0, 0, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s);
      else
        integral = (Xpa - (eta/p)*Xpq*s)*recurDKR(i-1, 0, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s)+ (1.0/(2.0*p))*(1-(eta/p)*s)*(i-1)*recurDKR(i-2, 0, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s);
    } else if (k == 0 && l == 0) {
      integral = (Xpb - (eta/p)*Xpq*s)*recurDKR(i, j-1, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s)+ (1.0/(2.0*p))*(1-(eta/p)*s)*(i*recurDKR(i-1, j-1, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (j-1)*recurDKR(i, j-2, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s));
    } else if (j== 0 && l == 0) {
      if (i > k) 
        integral = (Xpa - (eta/p)*Xpq*s)*recurDKR(i-1, 0, k, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*p))*(1-(eta/p)*s)*(i-1)*recurDKR(i-2, 0, k, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (s*k/(2.0*(p+q))) * recurDKR(i-1, 0, k-1, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s);
      else
        integral = (Xqc + (eta/q)*Xpq*s) * recurDKR(i, 0, k-1, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*(k-1)*recurDKR(i, 0, k-2, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (s*i/(2.0*(p+q))) * recurDKR(i-1, 0, k-1, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s);
    } else if (i == 0 && j == 0) {
      if (l == 0)  
        integral = (Xqc + (eta/q)*Xpq*s) * recurDKR(0, 0, k-1, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*(k-1)*recurDKR(0, 0, k-2, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s);
      else {
        if (k>l) 
        integral = (Xqc + (eta/q)*Xpq*s) * recurDKR(0, 0, k-1, l, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*((k-1)*recurDKR(0, 0, k-2, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + l*recurDKR(0, 0, k-1, l-1, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s));
        else
          integral = (Xqd + (eta/q)*Xpq*s) * recurDKR(0, 0, k, l-1, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*(k*recurDKR(0, 0, k-1, l-1, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (l-1)*recurDKR(0, 0, k, l-2, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s));
      } 
    } else if (l == 0) {
      if (k == 1)
        integral = (Xqc + (eta/q)*Xpq*s) * recurDKR(i, j, 0, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s)  + (s/(2.0*(p+q))) * (i*recurDKR(i-1,j,0,0,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s) + j*recurDKR(i,j-1,0,0,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s));
      else 
        integral = (Xqc + (eta/q)*Xpq*s) * recurDKR(i, j, k-1, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*(k-1)*recurDKR(i, j, k-2, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (s/(2.0*(p+q))) * (i*recurDKR(i-1,j,k-1,0,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s) + j*recurDKR(i,j-1,k-1,0,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s));
    } else {
      if (l == 1)  
        integral = (Xqd + (eta/q)*Xpq*s) * recurDKR(i, j, k, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*k*recurDKR(i, j, k-1, 0, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (s/(2.0*(p+q))) * (i*recurDKR(i-1,j,k,0,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s) + j*recurDKR(i,j-1,k,0,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s));
        else  
        integral = (Xqd + (eta/q)*Xpq*s) * recurDKR(i, j, k, l-1, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (1.0/(2.0*q))*(1-(eta/q)*s)*(k*recurDKR(i, j, k-1, l-1, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s) + (l-1)*recurDKR(i, j, k, l-2, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, s)) + (s/(2.0*(p+q))) * (i*recurDKR(i-1,j,k,l-1,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s) + j*recurDKR(i,j-1,k,l-1,Xpa,Xpb,Xqc,Xqd,Xpq,eta,p,q,s));
    }
  }
  return integral;

}
*/
double DKR::integrationBielec(int ix, int iy, int iz, int jx, int jy, int jz, int kx, int ky, int kz, int lx, int ly, int lz, double p, double q, double eta, double Xpa, double Ypa, double Zpa, double Xpb, double Ypb, double Zpb, double Xqc, double Yqc, double Zqc, double Xqd, double Yqd, double Zqd, double Xpq, double Ypq, double Zpq, double Rpq) {

  double factorX, factorY, factorZ, bielecTot, preTot;
/*
  bielecTot1 = 0.0;
  for (int n=0; n<Ns; n++) {
    result = 1;
    for (int i=0; i<3; i++)
     result *= recurDKR(ix, jx, kx, lx, Xpa, Xpb, Xqc, Xqd, Xpq, eta, p, q, var.R[n]);
    bielecTot1 += var.wi[n] * result;

  }
  
  bielecTot1 *= 0.5;
*/

  bielecTot = 0.0;
  for (int g=0; g<Ns; g++) {
    preTot = 0.0; 
    for (int t=0; t<=ix+jx; t++) 
      for (int u=0; u<=iy+jy; u++)
        for (int v=0; v<=iz+jz; v++)
          for (int ta=0; ta<=kx+lx; ta++) 
            for (int ua=0; ua<=ky+ly; ua++)
              for (int va=0; va<=kz+lz; va++) {
                factorX = pow(-1,ta) * MD.McMurchieCoef(ix, jx, t, p, Xpa, Xpb) * MD.McMurchieCoef(kx, lx, ta, q, Xqc, Xqd) * DKRrec(t+ta, Xpq, eta, var.R[g]);
                factorY = pow(-1,ua) * MD.McMurchieCoef(iy, jy, u, p, Ypa, Ypb) * MD.McMurchieCoef(ky, ly, ua, q, Yqc, Yqd) * DKRrec(u+ua, Ypq, eta, var.R[g]);

                factorZ = pow(-1,va) * MD.McMurchieCoef(iz, jz, v, p, Zpa, Zpb) * MD.McMurchieCoef(kz, lz, va, q, Zqc, Zqd) * DKRrec(v+va, Zpq, eta, var.R[g]);

                preTot += factorX * factorY * factorZ;
              }
    bielecTot += var.wi[g]*preTot;
  }

  bielecTot *= 0.5;

	return bielecTot;

}
/*
double DKR::integrationMonoelec() {

  double factorX, factorY, factorZ, bielecTot1, result, bielecTot, preTot;

  result = M_PI/dt.p;

  bielecTot = 0.0;
  for (int g=0; g<Ns; g++) {
    preTot = 0.0;
    for (int t=0; t<=dt.prim[0].ang[0]+dt.prim[1].ang[0]; t++)
      for (int u=0; u<=dt.prim[0].ang[1]+dt.prim[1].ang[1]; u++)
        for (int v=0; v<=dt.prim[0].ang[2]+dt.prim[1].ang[2]; v++) {
          factorX = McMurchieCoef(dt.prim[0].ang[0], dt.prim[1].ang[0], t, dt.p, dt.PA[0], dt.PB[0]) * DKRrec(t, dt.PQ[0], dt.p, var.R[g]);
          factorY = McMurchieCoef(dt.prim[0].ang[1], dt.prim[1].ang[1], u, dt.p, dt.PA[1], dt.PB[1]) * DKRrec(u, dt.PQ[1], dt.p, var.R[g]);
          factorZ = McMurchieCoef(dt.prim[0].ang[2], dt.prim[1].ang[2], v, dt.p, dt.PA[2], dt.PB[2]) * DKRrec(v, dt.PQ[2], dt.p, var.R[g]);

          preTot += factorX * factorY * factorZ;
          }

          bielecTot += var.wi[g]*preTot;
      }

  bielecTot *= result;

  cout <<"FINAL: " << setprecision(12) << bielecTot << endl;

return 0;

}
*/

