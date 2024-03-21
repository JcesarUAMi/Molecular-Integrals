#include <iostream>
#include <vector>
//#include <CL/sycl.hpp>
#include <cmath>
#include <fstream>
#include "readDataFiles.h"
using namespace std;

void readDataInput::coefficientsFinal(fstream& basisFile, fstream& xyzFile, fstream& movecsFile) {

  molec.readCoordinates(xyzFile);
  molec.getNumberOfDifAtoms();
  lastBasis.resize(molec.nNuc);
  basis.readBasisFile(basisFile, molec.nDifAtoms);
  coeff.readMovecsFile(movecsFile);
}

void readDataInput::nucleusBasis() {
  
  int atom;
  int j;
  numberOfMolOrb = coeff.c.nw_nbf;

  for (size_t i=0; i<molec.nNuc; i++) {
    atom = molec.geom[i].nAtomic;
    j = 0;
    while (j<molec.nDifAtoms) {
      if (atom == basis.basisPrimitiveFunctions[j].atomType) {
        lastBasis[i].atomType = atom;
        lastBasis[i].shellSize.resize(basis.basisPrimitiveFunctions[j].shellSize.size());
        lastBasis[i].functionType.resize(basis.basisPrimitiveFunctions[j].primitiveTypes.size());
        lastBasis[i].angularMoment.resize(basis.basisPrimitiveFunctions[j].totalAngMom.size());
        lastBasis[i].basisCoefficients.resize(basis.basisPrimitiveFunctions[j].primitiveCoefficients.size());
        lastBasis[i].basisExponents.resize(basis.basisPrimitiveFunctions[j].primitiveExponents.size());
        copy(begin(molec.geom[i].nCoor), end(molec.geom[i].nCoor), begin(lastBasis[i].nucleusCoord));  
        copy(begin(basis.basisPrimitiveFunctions[j].shellSize), end(basis.basisPrimitiveFunctions[j].shellSize), begin(lastBasis[i].shellSize));
        copy(begin(basis.basisPrimitiveFunctions[j].primitiveExponents), end(basis.basisPrimitiveFunctions[j].primitiveExponents), begin(lastBasis[i].basisExponents));
        copy(begin(basis.basisPrimitiveFunctions[j].primitiveCoefficients), end(basis.basisPrimitiveFunctions[j].primitiveCoefficients), begin(lastBasis[i].basisCoefficients));
        copy(begin(basis.basisPrimitiveFunctions[j].primitiveTypes), end(basis.basisPrimitiveFunctions[j].primitiveTypes), begin(lastBasis[i].functionType));
        copy(begin(basis.basisPrimitiveFunctions[j].totalAngMom), end(basis.basisPrimitiveFunctions[j].totalAngMom), begin(lastBasis[i].angularMoment));
        j = molec.nDifAtoms;
      } else
        j++;
    }
  }

  int cont, h;
  cont = 0;
  numberOfPrimFunc = 0;
  for (int i=0; i<lastBasis.size(); i++) {
    numberOfPrimFunc += lastBasis[i].basisExponents.size();  
    lastBasis[i].lastBasisCoeff.resize(coeff.c.nw_nbf); 
    for(int u=0; u<lastBasis[i].lastBasisCoeff.size(); u++) {
      lastBasis[i].lastBasisCoeff[u].bC.resize(lastBasis[i].shellSize.size());
      h = 0;
      for (int k=0; k<lastBasis[i].shellSize.size(); k++) {
        lastBasis[i].lastBasisCoeff[u].bC[k].resize(lastBasis[i].shellSize[k]);
        for (int l=0; l<lastBasis[i].shellSize[k]; l++) {
          lastBasis[i].lastBasisCoeff[u].bC[k][l] = coeff.molecularOrbitalData[u].orbitalCoeff[cont+k] * lastBasis[i].basisCoefficients[h];
          h++;
        }
      }
    }
   cont += lastBasis[i].shellSize.size();
  }

  moi.resize(numberOfPrimFunc);
}

void readDataInput::densityMatrix() {
  
  int cont, cont1, orb, tot;

  getNumOccOrb();
  tot = numberOfPrimFunc * numberOfPrimFunc;
  dE.resize(numberOfOccOrb*tot);
	pF.resize(numberOfOccOrb);
	cont = 0;
  for (int i=0; i<numberOfOccOrb; i++) 
    	for (int k=0; k<numberOfPrimFunc; k++) {
      	for (int j=0; j<numberOfPrimFunc; j++) {
        	dE[i*tot + numberOfPrimFunc*k+j] = 0.0;
      	}
    	}

  cont = 0;
	for (int k=0; k<numberOfOccOrb; k++)
		for (int i=0; i<coeff.c.nw_nbf; i++) {
			pF[k].dE.resize(coeff.c.nw_nbf);
			pF[k].dE[i].resize(coeff.c.nw_nbf);
			for (int j=0; j<coeff.c.nw_nbf; j++) {
				pF[k].dE[i][j] = 2.0*coeff.molecularOrbitalData[k].orbitalCoeff[i] * coeff.molecularOrbitalData[k].orbitalCoeff[j];	
//				cout << k << "   " << i << "  " << j << "   " << coeff.molecularOrbitalData[k].orbitalCoeff[i] << "   " << coeff.molecularOrbitalData[k].orbitalCoeff[j] << endl;
			}	
		}
 
	cont1 = 0;
	for (int j=0; j<lastBasis.size(); j++) 
		for (int k=0; k<lastBasis[j].shellSize.size(); k++) 
  		for (int l=0; l<lastBasis[j].shellSize[k]; l++) {
    		cont = 0;
    		for (int h=0; h<lastBasis.size(); h++) 
      		for (int m=0; m<lastBasis[h].shellSize.size(); m++)
        		for (int n=0; n<lastBasis[h].shellSize[m]; n++) { 
							for (int i=0; i<numberOfOccOrb; i++)  
          			dE[i*tot+cont1*numberOfPrimFunc+cont] = 2.0*lastBasis[h].lastBasisCoeff[i].bC[m][n] * lastBasis[j].lastBasisCoeff[i].bC[k][l];    
								cont++;
           }
				cont1++;
    		}

}  

void readDataInput::getNumOccOrb () {

  int occ;
  numberOfOccOrb = 0;
  numberOfVirOrb = 0; 
  for(int m=0; m<numberOfMolOrb; m++) {
    occ = coeff.molecularOrbitalData[m].orbitalOccupancy;
    if (occ != 0) {
      numberOfOccOrb++;
      Occ.push_back(m);
    } else {
      numberOfVirOrb++;
      Vir.push_back(m);
    }
  }
}

double readDataInput::generateMoi(double x, double y, double z) {
  
  double rr;
  double distx, disty, distz;
  double facx, facy, facz;
  double expo, prim, mo, den;
  int ang[3];
  int angMom;
  int funcType;
  int k;
  k = 0;
  for (int i=0; i<lastBasis.size(); i++)
    for (int j=0; j<lastBasis[i].basisExponents.size(); j++) {
      distx = x - lastBasis[i].nucleusCoord[0];
      disty = y - lastBasis[i].nucleusCoord[1];
      distz = z - lastBasis[i].nucleusCoord[2];
      rr = distx*distx + disty*disty + distz*distz;
      expo = exp(-lastBasis[i].basisExponents[j]*rr);
      angMom = lastBasis[i].angularMoment[j];
      funcType = lastBasis[i].functionType[j];
      basis.angularMoment(ang, angMom, funcType);
      facx = pow(distx, ang[0]); 
      facy = pow(disty, ang[1]); 
      facz = pow(distz, ang[2]); 

      moi[k] = facx*facy*facz*expo;
      k++;
    }

  den = 0.0;
  for (int i=0; i<numberOfOccOrb; i++) {
    mo = 0.0;
    k = 0;
    for (int h=0; h<lastBasis.size(); h++)
      for (int m=0; m<lastBasis[h].shellSize.size(); m++)
        for (int n=0; n<lastBasis[h].shellSize[m]; n++) {
          prim = moi[k]*lastBasis[h].lastBasisCoeff[i].bC[m][n];
          mo += prim;
          k++;
        }
        den += 2.0 * mo * mo;
      }

  return den;

}
