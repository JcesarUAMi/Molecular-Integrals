#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <string>
#include <array>
#include <vector>
using namespace std;

struct Geometry {
  array<double, 3> nCoor;
  int nAtomic;
};

class Molecule {
  private:
    string atomType;
  public:
    int nNuc;
    int nDifAtoms;
    int atomicNumber;
    void getNumberOfDifAtoms();
    void readCoordinates(fstream&);
    void getAtomicNumber (string);
    vector<Geometry> geom;
    vector<int> numDifAtoms;

};


#endif
