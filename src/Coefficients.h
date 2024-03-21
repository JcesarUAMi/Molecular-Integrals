#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

struct coeffData {
  vector<double> orbitalCoeff;
  double orbitalOccupancy;
  double orbitalEnergie;
};

struct readDataMovecs {
    int ignoreInt;
    char ignoreChar[256];
    int titleLength;
    int basisTitleLength;
    int nw_sets;
    int nw_nbf;
    int nw_nmo;
    double *nw_en;
    double *nw_co;
    double *occ;
};

class Coefficients {
  public:
    ~Coefficients();
    vector<coeffData> molecularOrbitalData;
    void readMovecsFile(fstream&);
    readDataMovecs c;  

};

#endif
