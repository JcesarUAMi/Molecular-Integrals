#ifndef BASIS_H
#define BASIS_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "Geometry.h"
using namespace std;

struct basisPrimitiveData {
  int atomType;
  vector<int> primitiveTypes;
  vector<double> primitiveExponents;
  vector<double> primitiveCoefficients;
  vector<int> shellSize;
  vector<int> totalAngMom;
  vector<double> normConst;
};

class Basis : public Molecule {
  public:
    int SoC;
    string type = "Cartesian";
    vector<basisPrimitiveData> basisPrimitiveFunctions;
    void readBasisFile (fstream&, int);
    void angularMoment(int*, int, int);
    double factorialDouble (int);
  private:
    string gaussianType;
    double normalizeCartesian (double, int*);
    double normalizeSpherical (double, int);
    void normalize(int);
    int getSymmetry (string, int);
};

#endif
  


 
