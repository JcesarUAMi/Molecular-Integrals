#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
//#include <CL/sycl.hpp>
#include "McMurchie.h"
#include "BoysFunction.h"
using namespace std;

void McMD::LectureBoysValues (fstream& readFromFile, fstream& readFromFile1) {

  double ta, tb, tc, td;
  string file = "FinalCubicSplineData.dat";
  ifstream boysFile(file.c_str(), ios::binary | ios::in);
  if (!boysFile) {
    cerr << "File couldn't be opened" << endl;
    exit (EXIT_FAILURE);
  }

  for (int i=0; i<10; i++) {
    boysF.btv[i].a.resize(300001);
    boysF.btv[i].b.resize(300001);
    boysF.btv[i].c.resize(300001);
    boysF.btv[i].d.resize(300001);
    for(int j=0; j<300001; j++) {
      boysFile.read(reinterpret_cast<char *>(&ta), sizeof(ta));
      boysFile.read(reinterpret_cast<char *>(&tb), sizeof(tb));
      boysFile.read(reinterpret_cast<char *>(&tc), sizeof(tc));
      boysFile.read(reinterpret_cast<char *>(&td), sizeof(td));
      boysF.btv[i].a[j] = ta;
      boysF.btv[i].b[j] = tb;
      boysF.btv[i].c[j] = tc;
      boysF.btv[i].d[j] = td;
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
    boysF.btvNum[i].a.resize(300001);
    boysF.btvNum[i].b.resize(300001);
    boysF.btvNum[i].c.resize(300001);
    boysF.btvNum[i].d.resize(300001);
    for(int j=0; j<300001; j++) {
      boysFile1.read(reinterpret_cast<char *>(&ta), sizeof(ta));
      boysFile1.read(reinterpret_cast<char *>(&tb), sizeof(tb));
      boysFile1.read(reinterpret_cast<char *>(&tc), sizeof(tc));
      boysFile1.read(reinterpret_cast<char *>(&td), sizeof(td));
      boysF.btvNum[i].a[j] = ta;
      boysF.btvNum[i].b[j] = tb;
      boysF.btvNum[i].c[j] = tc;
      boysF.btvNum[i].d[j] = td;
    }
  }

  boysFile1.close();

}

double McMD::McMurchieIntegral (int n, double redExpo, double Rpq) {

  double integral;
  double factor;
  factor = redExpo * Rpq;
  if (n == 0) {
  	integral = boysF.boysZeroDegree<double>(factor);
	} else {
  	integral = pow(-2.0*redExpo, n); 
    integral *= boysF.BoysFunction<double>(n, factor);
		}

  return integral;

}

double McMD::McMurchieCoefOC(int la, int lb, int t, double p) {

    if (t < 0 || t > la+lb)
    return 0.0;
  if (t == 0) {
    if (la == 0 && lb == 0)
      return 1.0;
    else
      if (la > lb)
        return McMurchieCoefOC(la-1, lb, 1, p);
      else
        return McMurchieCoefOC(la, lb-1, 1, p);
  } else
    if (t == 1)
      return (1 / (2 * p)) * (la * McMurchieCoefOC(la-1, lb, 0, p) + lb * McMurchieCoefOC(la, lb-1, 0, p));
    else
      return (1 / (2 * p * t)) * (la * McMurchieCoefOC(la-1, lb, t-1, p) + lb * McMurchieCoefOC(la, lb-1, t-1, p));

}

double McMD::McMurchieCoef (int la, int lb, int t, double p, double Xpa, double Xpb) {
  
  if (t < 0 || t > la+lb)
    return 0.0;
  if (t == 0) {
    if (la == 0 && lb == 0)
      return 1.0;
    else 
      if (la > lb)
        return Xpa * McMurchieCoef(la-1, lb, 0, p, Xpa, Xpb) + McMurchieCoef(la-1, lb, 1, p, Xpa, Xpb);
      else
        return Xpb * McMurchieCoef(la, lb-1, 0, p, Xpa, Xpb) + McMurchieCoef(la, lb-1, 1, p, Xpa, Xpb);
  } else 
    if (t == 1)
      return (1 / (2 * p)) * (la * McMurchieCoef(la-1, lb, 0, p, Xpa, Xpb) + lb * McMurchieCoef(la, lb-1, 0, p, Xpa, Xpb));
    else
      return (1 / (2 * p * t)) * (la * McMurchieCoef(la-1, lb, t-1, p, Xpa, Xpb) + lb * McMurchieCoef(la, lb-1, t-1, p, Xpa, Xpb));

}

double McMD::McMurchieRecursionOC (int t, int u, int v, int n, double redExpo) {

  double integral;
  if (t < 0 || u < 0 || v < 0) {
    integral = 0.0;
  } else {
    if (t == 0 && u == 0 && v == 0) {
      if(n<0) n = 0;
      integral = McMurchieIntegral(n, redExpo, 0.0);
    } else {
      if (t > 0)
        integral = (t-1)*McMurchieRecursionOC(t-2, u, v, n-1, redExpo);
      else
        if (u > 0)
          integral = (u-1)*McMurchieRecursionOC(t, u-2, v, n-1, redExpo);
        else
          if (v > 0)
            integral = (v-1)*McMurchieRecursionOC(t, u, v-2, n-1, redExpo);
    }
  }

  return integral;
}

double McMD::McMurchieRecursion (int t, int u, int v, int n, double Rpq, double Xpq, double Ypq, double Zpq, double redExpo) {

  double integral;
 
 	if (t < 0 || u < 0 || v < 0) {
    integral = 0.0;
 	} else {
    if (t == 0 && u == 0 && v == 0) {
			if(n<0) n = 0;
      integral = McMurchieIntegral(n, redExpo, Rpq);
  	} else {
      if (t > 0)
        integral = ((t-1)*McMurchieRecursion(t-2, u, v, n-1, Rpq, Xpq, Ypq, Zpq, redExpo) + Xpq*McMurchieRecursion(t-1, u, v, n, Rpq, Xpq, Ypq, Zpq, redExpo));
      else
        if (u > 0)
          integral = ((u-1)*McMurchieRecursion(t, u-2, v, n-1, Rpq, Xpq, Ypq, Zpq, redExpo) + Ypq*McMurchieRecursion(t, u-1, v, n, Rpq, Xpq, Ypq, Zpq, redExpo));
        else
          if (v > 0) 
            integral = (v-1)*McMurchieRecursion(t, u, v-2, n-1, Rpq, Xpq, Ypq, Zpq, redExpo) + Zpq*McMurchieRecursion(t, u, v-1, n, Rpq, Xpq, Ypq, Zpq, redExpo);
					
    }

  }
  return integral;
}

double McMD::ObaraIntegral (int n, double redExpo, double Rpq) {

  double integral;
  double factor;
	
  factor = redExpo * Rpq;

  if (n == 0) 
    integral = boysF.boysZeroDegree<double>(factor);
  else 
    integral = boysF.BoysFunction<double>(n, factor);
	
  return integral;

}

double McMD::newObara (int ix, int iy, int iz, int jx, int jy, int jz, int kx, int ky, int kz, int lx, int ly, int lz, int n, double p, double q, double eta, double Xpa, double Ypa, double Zpa, double Xpb, double Ypb, double Zpb, double Xqc, double Yqc, double Zqc, double Xqd, double Yqd, double Zqd, double Xpq, double Ypq, double Zpq, double Rpq) {

  double integral;
  int tx, ty, tz;

  tx = ix + jx + kx + lx;
  ty = iy + jy + ky + ly;
  tz = iz + jz + kz + lz;

  if (tx == 0 && ty == 0 && tz == 0)
    integral = ObaraIntegral(n, eta, Rpq);
   else if (tx > 0) {

    if (kx == 0 && lx == 0) {
      if (jx == 0) {
        if (ix == 1)
          integral = Xpa * newObara(0, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p) * Xpq * newObara(0, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      else if (ix > 1) {
        n--;
        integral = Xpa*newObara(ix-1, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p)*Xpq*newObara(ix-1, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (ix-1)/(2.0*p) * (newObara(ix-2, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p) * newObara(ix-2, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else {
      integral = newObara(ix+1, iy, iz, jx-1, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Xpb - Xpa)*newObara(ix, iy, iz, jx-1, jy, jz, 0, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      }
    } else if (ix == 0 && jx == 0) {
    if (lx == 0) {
      if (kx == 1)
        integral = Xqc * newObara(0, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Xpq * newObara(0, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      else if (kx > 1) {
        n--;
        integral = Xqc*newObara(0, iy, iz, 0, jy, jz, kx-1, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q)*Xpq*newObara(0, iy, iz, 0, jy, jz, kx-1, ky, kz, 0, ly, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (kx-1)/(2.0*q) * (newObara(0, iy, iz, 0, jy, jz, kx-2, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/q) * newObara(0, iy, iz, 0, jy, jz, kx-2, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else {
      integral = newObara(0, iy, iz, 0, jy, jz, kx+1, ky, kz, lx-1, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Xqd - Xqc)*newObara(0, iy, iz, 0, jy, jz, kx, ky, kz, lx-1, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      }
    } else if (lx == 0 && jx == 0) {
      n--;
      if (kx > 0) {
        if (kx == 1) {
          integral =  Xqc * newObara(ix, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Xpq * newObara(ix, iy, iz, jx, jy, jz, 0, ky, kz, 0, ly, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (ix/(2.0*(p+q))) * newObara(ix-1, iy, iz, 0, jy, jz, 0, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
        } else
          integral =  Xqc * newObara(ix, iy, iz, 0, jy, jz, kx-1, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Xpq * newObara(ix, iy, iz, 0, jy, jz, kx-1, ky, kz, 0, ly, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (ix/(2.0*(p+q))) * newObara(ix-1, iy, iz, 0, jy, jz, kx-1, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + ((kx-1)/(2.0*q)) * (newObara(ix, iy, iz, 0, jy, jz, kx-2, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/q) * newObara(ix, iy, iz, 0, jy, jz, kx-2, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else if (lx == 0)
      integral = newObara(ix+1, iy, iz, jx-1, jy, jz, kx, ky, kz, 0, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Xpb - Xpa)*newObara(ix, iy, iz, jx-1, jy, jz, kx, ky, kz, 0, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
    else
      integral = newObara(ix, iy, iz, jx, jy, jz, kx+1, ky, kz, lx-1, ly, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Xqd - Xqc)*newObara(ix, iy, iz, jx, jy, jz, kx, ky, kz, lx-1, ly, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);

  } else if (ty > 0) {
    if (ky == 0 && ly == 0) {
      if (jy == 0) {
        if (iy == 1)
          integral = Ypa * newObara(ix, 0, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p) * Ypq * newObara(ix, 0, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      else if (iy > 1) {
        n--;
        integral = Ypa*newObara(ix, iy-1, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p)*Ypq*newObara(ix, iy-1, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (iy-1)/(2.0*p) * (newObara(ix, iy-2, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p) * newObara(ix, iy-2, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else {
      integral = newObara(ix, iy+1, iz, jx, jy-1, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Ypb - Ypa)*newObara(ix, iy, iz, jx, jy-1, jz, kx, 0, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      }
    } else if (iy == 0 && jy == 0) {
    if (ly == 0) {
      if (ky == 1)
        integral = Yqc * newObara(ix, 0, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Ypq * newObara(ix, 0, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      else if (ky > 1) {
        n--;
        integral = Yqc*newObara(ix, 0, iz, jx, 0, jz, kx, ky-1, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q)*Ypq*newObara(ix, 0, iz, jx, 0, jz, kx, ky-1, kz, lx, 0, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (ky-1)/(2.0*q) * (newObara(ix, 0, iz, jx, 0, jz, kx, ky-2, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/q) * newObara(ix, 0, iz, jx, 0, jz, kx, ky-2, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else {
      integral = newObara(ix, 0, iz, jx, 0, jz, kx, ky+1, kz, lx, ly-1, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Yqd - Yqc)*newObara(ix, 0, iz, jx, 0, jz, kx, ky, kz, lx, ly-1, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      }
    } else if (ly == 0 && jy == 0) {
      n--;
      if (ky > 0) {
        if (ky == 1) {
          integral =  Yqc * newObara(ix, iy, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Ypq * newObara(ix, iy, iz, jx, jy, jz, kx, 0, kz, lx, 0, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (iy/(2.0*(p+q))) * newObara(ix, iy-1, iz, jx, 0, jz, kx, 0, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
        } else
          integral =  Yqc * newObara(ix, iy, iz, jx, 0, jz, kx, ky-1, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Ypq * newObara(ix, iy, iz, jx, 0, jz, kx, ky-1, kz, lx, 0, lz, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (iy/(2.0*(p+q))) * newObara(ix, iy-1, iz, jx, 0, jz, kx, ky-1, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + ((ky-1)/(2.0*q)) * (newObara(ix, iy, iz, jx, 0, jz, kx, ky-2, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/q) * newObara(ix, iy, iz, jx, 0, jz, kx, ky-2, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else if (ly == 0)
      integral = newObara(ix, iy+1, iz, jx, jy-1, jz, kx, ky, kz, lx, 0, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Ypb - Ypa)*newObara(ix, iy, iz, jx, jy-1, jz, kx, ky, kz, lx, 0, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
    else
      integral = newObara(ix, iy, iz, jx, jy, jz, kx, ky+1, kz, lx, ly-1, lz, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Yqd - Yqc)*newObara(ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly-1, lz, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);

  } else if (tz > 0) {

   if (kz == 0 && lz == 0) {
      if (jz == 0) {
        if (iz == 1)
          integral = Zpa * newObara(ix, iy, 0, jx, jy, 0, kx, ky, 0, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p) * Zpq * newObara(ix, iy, 0, jx, jy, 0, kx, jy, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      else if (iz > 1) {
        n--;
        integral = Zpa*newObara(ix, iy, iz-1, jx, jy, 0, kx, ky, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p)*Zpq*newObara(ix, iy, iz-1, jx, jy, 0, kx, ky, 0, lx, ly, 0, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (iz-1)/(2.0*p) * (newObara(ix, iy, iz-2, jx, jy, 0, kx, ky, 0, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/p) * newObara(ix, iy, iz-2, jx, jy, 0, kx, ky, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else {
      integral = newObara(ix, iy, iz+1, jx, jy, jz-1, kx, ky, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Zpb - Zpa)*newObara(ix, iy, iz, jx, jy, jz-1, kx, ky, 0, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      }
    } else if (iz == 0 && jz == 0) {
    if (lz == 0) {
      if (kz == 1)
        integral = Zqc * newObara(ix, iy, 0, jx, jy, 0, kx, ky, 0, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Zpq * newObara(ix, iy, 0, jx, jy, 0, kx, ky, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      else if (kz > 1) {
        n--;
        integral = Zqc*newObara(ix, iy, 0, jx, jy, 0, kx, ky, kz-1, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q)*Zpq*newObara(ix, iy, 0, jx, jy, 0, kx, ky, kz-1, lx, ly, 0, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (kz-1)/(2.0*q) * (newObara(ix, iy, 0, jx, jy, 0, kx, ky, kz-2, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/q) * newObara(ix, iy, 0, jx, jy, 0, kx, ky, kz-2, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else {
      integral = newObara(ix, iy, 0, jx, jy, 0, kx, ky, kz+1, lx, ly, lz-1, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Zqd - Zqc)*newObara(ix, iy, 0, jx, jy, 0, kx, ky, kz, lx, ly, lz-1, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
      }
    } else if (lz == 0 && jz == 0) {
      n--;
      if (kz > 0) {
        if (kz == 1) {
          integral =  Zqc * newObara(ix, iy, iz, jx, jy, 0, kx, ky, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Zpq * newObara(ix, iy, iz, jx, jy, jz, kx, ky, 0, lx, ly, 0, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (iz/(2.0*(p+q))) * newObara(ix, iy, iz-1, jx, jy, 0, kx, ky, 0, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
        } else
          integral =  Zqc * newObara(ix, iy, iz, jx, jy, 0, kx, ky, kz-1, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (eta/q) * Zpq * newObara(ix, iy, iz, jx, jy, 0, kx, ky, kz-1, lx, ly, 0, n+1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (iz/(2.0*(p+q))) * newObara(ix, iy, iz-1, jx, jy, 0, kx, ky, kz-1, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + ((kz-1)/(2.0*q)) * (newObara(ix, iy, iz, jx, jy, 0, kx, ky, kz-2, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) - (eta/q) * newObara(ix, iy, iz, jx, jy, 0, kx, ky, kz-2, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq));
      }
    } else if (lz == 0)
      integral = newObara(ix, iy, iz+1, jx, jy, jz-1, kx, ky, kz, lx, ly, 0, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Zpb - Zpa)*newObara(ix, iy, iz, jx, jy, jz-1, kx, ky, kz, lx, ly, 0, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);
    else
      integral = newObara(ix, iy, iz, jx, jy, jz, kx, ky, kz+1, lx, ly, lz-1, n, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq) + (Zqd - Zqc)*newObara(ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly, lz-1, n-1, p, q, eta, Xpa, Ypa, Zpa, Xpb, Ypb, Zpb, Xqc, Yqc, Zqc, Xqd, Yqd, Zqd, Xpq, Ypq, Zpq, Rpq);


  }

  return integral;

}

double McMD::TraslapeNumericalIntegral(int la, int lb, double Xpa, double Xpb, double q) {
    int i;
    double integral = 0.0, sq;
    double x, xk[12], wk[12];
    sq = sqrt(q);

  xk[0]  = -3.889724897869781919272;     xk[1]  = -3.020637025120889771711;
  xk[2]  = -2.279507080501059900188;     xk[3]  = -1.59768263515260479671;
  xk[4]  = -0.9477883912401637437046;    xk[5]  = -0.314240376254359111277;
  xk[6]  =  3.889724897869781919272;     xk[7]  =  3.020637025120889771711;
  xk[8]  =  2.279507080501059900188;     xk[9]  =  1.59768263515260479671;
  xk[10] = 0.9477883912401637437046;     xk[11] =  0.314240376254359111277;
  wk[0] = 2.65855168435630160602E-7;     wk[1]  = 8.5736870435878586546E-5;
  wk[2] = 0.00390539058462906185999;     wk[3]  = 0.05160798561588392999187;
  wk[4] = 0.2604923102641611292334;      wk[5]  = 0.5701352362624795783471;
  wk[6] = 2.65855168435630160602E-7;     wk[7]  = 8.5736870435878586546E-5;
  wk[8] = 0.00390539058462906185999;     wk[9]  = 0.05160798561588392999187;
  wk[10] = 0.2604923102641611292334;     wk[11] = 0.5701352362624795783471;

  if (la == 0 && lb == 0)
    return sqrt(M_PI) / sq;
  else {
    for (i = 0; i < 12; i++) {
      x = xk[i]/sq;
      integral += wk[i] * pow(x + Xpa, la) * pow(x + Xpb, lb);
    }
  return integral/sq;
  }

}

double McMD::overlapMcMurchie(int la, int lb, double expoTot, double Xpa, double Xpb) {

  double integral;

  if (la == 0 && lb == 0)
    integral = sqrt(M_PI/expoTot);
  else
    integral = sqrt(M_PI/expoTot) * McMurchieCoef(la, lb, 0, expoTot, Xpa, Xpb);

  return integral;

}




