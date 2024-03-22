#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <CL/sycl.hpp>
using namespace std;

#include "Coefficients.h"

void Coefficients::readMovecsFile(fstream& readFromFile) {

//  string file = "h2dif.movecs";
  string file = "decanol.movecs";
//  string file = "../h2oaug.movecs";
//  string file = "heaug.movecs";
//  string file = "he.movecs";
//  string file = "coohdif.movecs";
//  string file = "hepvt.movecs";
//  string file = "h2osto.movecs";
//  string file = "h2o631.movecs";
//  string file = "o2aug.movecs";
//  string file = "benzene.movecs";
  ifstream movecs(file.c_str(), ios::in | ios::binary);

  if (!movecs) {
    cerr << "The Coefficients movecs file can't be open." << "\n";
    exit(EXIT_FAILURE);
  }

  // get the calculation info.
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read(c.ignoreChar, c.ignoreInt);
  c.ignoreChar[c.ignoreInt] = '\0';
  movecs.read((char*)&c.ignoreInt, 4);

  //get calculation type
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read(c.ignoreChar, c.ignoreInt);
  c.ignoreChar[c.ignoreInt] = '\0';
  movecs.read((char*)&c.ignoreInt, 4);

  //get title length
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)&c.titleLength, c.ignoreInt);
  movecs.read((char*)&c.ignoreInt, 4);

  //title
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read(c.ignoreChar, c.ignoreInt);
  c.ignoreChar[c.ignoreInt] = '\0';
  movecs.read((char*)&c.ignoreInt, 4);

  //basis name length
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)&c.basisTitleLength, c.ignoreInt);
  movecs.read((char*)&c.ignoreInt, 4);

  //basis name
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read(c.ignoreChar, c.ignoreInt);
  c.ignoreChar[c.ignoreInt] = '\0';
  movecs.read((char*)&c.ignoreInt, 4);

  //nw_sets
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)&c.nw_sets, c.ignoreInt);
  movecs.read((char*)&c.ignoreInt, 4);

  //nw_nbf
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)&c.nw_nbf, c.ignoreInt);
  movecs.read((char*)&c.ignoreInt, 4);

  //nw_nmo
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)&c.nw_nmo, c.ignoreInt);
  movecs.read((char*)&c.ignoreInt, 4);
  
  cout.flush();

  molecularOrbitalData.resize(c.nw_nbf);
  c.occ = new double[c.nw_nbf];
  c.nw_en = new double[c.nw_nbf];
  c.nw_co = new double[c.nw_nbf*c.nw_nbf];

  //occupancy
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)c.occ, c.ignoreInt);
  
  for(int i=0;i<c.nw_nbf;i++) { 
//    cout << i << " oillaaa: " << c.occ[i] << endl;
    molecularOrbitalData[i].orbitalOccupancy = c.occ[i];
  }
  
  movecs.read((char*)&c.ignoreInt, 4);

  //energies
  movecs.read((char*)&c.ignoreInt, 4);
  movecs.read((char*)c.nw_en, c.ignoreInt);
  
  for(int i=0;i<c.nw_nbf;i++) {  
    molecularOrbitalData[i].orbitalEnergie = c.nw_en[i];
//    cout << i << " baiiii: " << c.nw_en[i] << endl;
  }
  
  movecs.read((char*)&c.ignoreInt, 4);

  for (int i=0; i<c.nw_nbf; i++) {
    movecs.read((char*)&c.ignoreInt, 4);
    movecs.read((char*)c.nw_co, c.ignoreInt);
    molecularOrbitalData[i].orbitalCoeff.resize(c.nw_nbf);
    for (int j=0; j<c.nw_nbf; j++) 
      molecularOrbitalData[i].orbitalCoeff[j] = c.nw_co[j];
    
    movecs.read((char*)&c.ignoreInt, 4);
  }

  movecs.close();

/*
  for (int i=0; i<c.nw_nbf; i++) {
    cout << "FOR ORBITAL: " << i << endl;
    cout << "Occupancy: " << molecularOrbitalData[i].orbitalOccupancy << endl;
    cout << "Energie: " << molecularOrbitalData[i].orbitalEnergie << endl;
    cout << "Orbital Coeffs: " << endl;
    for (int j=0; j<c.nw_nbf; j++) 
      cout << j << "   " << molecularOrbitalData[i].orbitalCoeff[j] << endl;
  }
*/
}

Coefficients::~Coefficients() {
  delete[] c.nw_en;
  delete[] c.nw_co;
  delete[] c.occ;

}

