#ifndef MCMD_H
#define MCMD_H

#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "BoysFunction.h"
using namespace std;

class McMD {
  private:
    Boys<double> boysF;
  public:
    double McMurchieRecursion(int, int, int, int, double, double, double, double, double);
    double McMurchieCoef(int, int, int, double, double, double);
    double McMurchieRecursionOC(int, int, int, int, double);
    double McMurchieCoefOC(int, int, int, double);
		void LectureBoysValues(fstream&, fstream&);
		double ObaraRecursion(int, int, int, int, int, double, double, double, double, double, double, double, double, double);
		double TraslapeNumericalIntegral (int, int, double, double, double);
    double overlapMcMurchie(int, int, double, double, double);
		double newObara (int, int, int, int, int, int, int, int, int, int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
    double ObaraIntegral(int, double, double);
    double McMurchieIntegral(int, double, double);
};

#endif
