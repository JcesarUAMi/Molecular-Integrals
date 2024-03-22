#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <algorithm>
//#include <CL/sycl.hpp>
using namespace std;
#include "Geometry.h"

void Molecule::readCoordinates(fstream& fileXYZ) {
  
  double coor[3];
  const double convertAngBohr = 1.889725989;
//  string file = "h2oSto.xyz";
//  string file = "h2dif.xyz";
//  string file = "h2odif.xyz";
//  string file = "O2.xyz";
//  string file = "h2o631.xyz";
  string file = "decanol.xyz";
//  string file = "cooh.xyz";
//  string file = "h2oaug.xyz";
//  string file = "benzene.xyz";
  ifstream geometry(file.c_str(), ios::in);
  
  if (!geometry) {
    cerr << "The geometry xyz file can't be open." << endl;
    exit(EXIT_FAILURE);
  }
  
  geometry >> nNuc;
  geometry.ignore();
  geom.resize(nNuc);
//  cout << "there are " << nNuc << " atoms in total." << endl;
  for(int i=0; i<nNuc; ++i) {
    coor[0] = coor[1] = coor[2] = 0.0;
    atomicNumber = 0;
    geometry >> atomType >> coor[0] >> coor[1] >> coor[2];
    getAtomicNumber(atomType);
    geom[i].nAtomic = atomicNumber;
    geom[i].nCoor[0] = coor[0] * convertAngBohr;
    geom[i].nCoor[1] = coor[1] * convertAngBohr;
    geom[i].nCoor[2] = coor[2] * convertAngBohr;
  }
  geometry.close();
/*  
  for (int j=0; j <nNuc; ++j) {
    cout << "the atomic Number is " << geom[j].nAtomic
    << " with coordinates: "
    << geom[j].nCoor[0] << '\t' << geom[j].nCoor[1] << '\t' << geom[j].nCoor[2] << endl;
  }
*/  
}

void Molecule::getAtomicNumber (string atomType) {
  array<string, 118> tableElements = {
"H",                                                                                 "He",
"Li","Be",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
"Na","Mg",                                                  "Al","Si","P" ,"S" ,"Cl","Ar",
"K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
"Rb","Sr","Y" ,"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe",
"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
"Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
"Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf",
"Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
 };
  auto location = find(tableElements.cbegin(), tableElements.cend(), atomType);
  if (location != tableElements.cend())
    atomicNumber = location - tableElements.cbegin() + 1;
  else
    cout << "\nCheck the atom name in the xyz file." << endl;
}

void Molecule::getNumberOfDifAtoms () {
  
  int t1, t2, j;
 
  numDifAtoms.resize(1);
  numDifAtoms[0] = geom[0].nAtomic;
	j = 0;
	if (nNuc > 1) {
  	for (int i=0; i<nNuc; i++) {
    	t1 = geom[i].nAtomic; 
    	j = 0;
    	while (j < numDifAtoms.size()) {
      	t2 = numDifAtoms[j];
      	if (t1 == t2) 
        	j = nNuc;
    	 	else 
        	j++;      
    	}
    if (j == numDifAtoms.size())
      numDifAtoms.push_back(t1);
  	}
	}

  nDifAtoms = numDifAtoms.size();
}
