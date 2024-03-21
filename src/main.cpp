#include <iostream>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <CL/sycl.hpp>
#include <bits/stdc++.h>
#include "Coefficients.h"
#include "Geometry.h"
#include "Basis.h"
#include "readDataFiles.h"
#include "MatElements.h"
#include "BoysFunction.h"
#include "McMurchie.h"
using namespace std;

int main () {

  fstream xyz;
  fstream movecs;
  fstream basis;
  fstream boys;
  fstream boys1;
	MatElements ME;
	ME.Data(basis, xyz, movecs, boys, boys1);

//  ME.detectDevice();
  auto begin = chrono::high_resolution_clock::now();
  ME.matrixElements(boys, boys1);
  auto end = chrono::high_resolution_clock::now();
  auto elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

  cout << "Time measured to solve the matrix elements: " <<  elapsed.count() * 1e-9 << "\n";
 
  begin = chrono::high_resolution_clock::now();
  ME.matrixBuilding();
  end = chrono::high_resolution_clock::now();
  elapsed = chrono::duration_cast<std::chrono::nanoseconds>(end - begin);

  cout << "Time measured to build the Fock matrix: " << elapsed.count() * 1e-9 << "\n";
 
  return 0;

}

