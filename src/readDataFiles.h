#ifndef READ_H
#define READ_H

#include "Basis.h"
#include "Geometry.h"
#include "Coefficients.h"
#include <iostream>
#include <vector>
#include <array>
#include <iomanip>
#include <cstdlib>
using namespace std;

struct basisCoeff {
  vector<vector<double>> bC;
};

struct finalBasisData { 
  vector<double> basisExponents;
  vector<double> basisCoefficients;
  vector<basisCoeff> lastBasisCoeff;
  vector<int> shellSize;
  vector<int> functionType;
  vector<int> angularMoment;
  array<double, 3> nucleusCoord;
  int atomType;
	double charge;
}; 

struct densityElements {
  vector<vector<double>> dE;
};

class readDataInput {
  public:
    int numberOfMolOrb;
    int numberOfOccOrb;
    int numberOfVirOrb;
    int numberOfPrimFunc;
    vector<double> moi;
    vector<double> dE;
		vector<densityElements> pF;
    vector<finalBasisData> lastBasis;
    double generateMoi(double, double, double);
    void densityMatrix();
    void coefficientsFinal(fstream&, fstream&, fstream&);
    void nucleusBasis();
    void getNumOccOrb();
    vector<int> Occ;
    vector<int> Vir;
  private:
    Basis basis;
    Molecule molec;
    Coefficients coeff;

};

#endif


