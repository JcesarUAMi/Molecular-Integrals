#include <iostream>
#include <CL/sycl.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <iomanip>
#include "MatElements.h"
using namespace std;
using namespace sycl;


void MatElements::OrbitalData(int num) {
  
  int atomCenter2, atomCenter1, angMom1, angMom2, ang1[3], ang2[3], angs[3], angT;
  double expo1, expo2, Ax, Ay, Az, Bx, By, Bz, expoTot, mu, Px, Py, Pz, AB, Kab, PA[3], PB[3];
  int u, la, ni, lb;
  int n = num * (num - 1 ) / 2;
  ni = num * (num + 1 ) / 2;
  oEi.resize(n);
  oETS.resize(num);
  schF.resize(num*num);
  u = 0;
  la = 0;
  for (int i=0; i<num; i++) {
    expo1 = aO[i].expo;
    Ax = aO[i].coord[0];
    Ay = aO[i].coord[1];
    Az = aO[i].coord[2];
    atomCenter1 = aO[i].atomCenter;
    angMom1 = aO[i].angMom;
    for (int ip=0; ip<3; ip++)
      ang1[ip] = aO[i].ang[ip];

    for (int j=i; j<num; j++) {
      expo2 = aO[j].expo;
      Bx = aO[j].coord[0];
      By = aO[j].coord[1];
      Bz = aO[j].coord[2];
      atomCenter2 = aO[j].atomCenter;
      angMom2 = aO[j].angMom;
      for (int ip=0; ip<3; ip++)
        ang2[ip] = aO[j].ang[ip];

      expoTot = expo1 + expo2;
      angT = angMom1 + angMom2;
      angs[0] = ang1[0] + ang2[0];
      angs[1] = ang1[1] + ang2[1];
      angs[2] = ang1[2] + ang2[2];


      mu = (expo1 * expo2) / expoTot;
      Px = (expo1 * Ax + expo2 * Bx) / expoTot;
      Py = (expo1 * Ay + expo2 * By) / expoTot;
      Pz = (expo1 * Az + expo2 * Bz) / expoTot;

      AB = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
      Kab = -mu * AB;
  
      PA[0] = Px - Ax;
      PA[1] = Py - Ay;
      PA[2] = Pz - Az;

      PB[0] = Px - Bx;
      PB[1] = Py - By;
      PB[2] = Pz - Bz;

      la = num*i + j;
      lb = num*j + i;

      if (i != j ) {
        oEi[u].ind_a = i;
        oEi[u].ind_b = j;
        u++;
      } else {
        oETS[i].ind_a = i;
        oETS[i].ind_b = j;
      }

      schF[la].ind_a = schF[lb].ind_a = i;
      schF[la].ind_b = schF[lb].ind_b = j;
      schF[la].atmCenter1 = schF[lb].atmCenter1 = atomCenter1;
      schF[la].atmCenter2 = schF[lb].atmCenter2 = atomCenter2;
      schF[la].expo1 = schF[lb].expo1 = expo1;
      schF[la].expo2 = schF[lb].expo2 = expo2;
      schF[la].expoTot = schF[lb].expoTot = expoTot;
      schF[la].ext = schF[lb].ext = sqrt(2.0/(expoTot))*icef;
      schF[la].mu = schF[lb].mu = mu;
      schF[la].P[0] = schF[lb].P[0] = Px;
      schF[la].P[1] = schF[lb].P[1] = Py;
      schF[la].P[2] = schF[lb].P[2] = Pz;
      schF[la].angT = schF[lb].angT = angT;

      schF[la].angs[0] = schF[lb].angs[0] = angs[0];
      schF[la].angs[1] = schF[lb].angs[1] = angs[1];
      schF[la].angs[2] = schF[lb].angs[2] = angs[2];

      schF[la].ang1[0] = schF[lb].ang1[0] = ang1[0];
      schF[la].ang1[1] = schF[lb].ang1[1] = ang1[1];
      schF[la].ang1[2] = schF[lb].ang1[2] = ang1[2];

      schF[la].ang2[0] = schF[lb].ang2[0] = ang2[0];
      schF[la].ang2[1] = schF[lb].ang2[1] = ang2[1];
      schF[la].ang2[2] = schF[lb].ang2[2] = ang2[2];

      schF[la].AB = schF[lb].AB = AB;
      schF[la].Kab = schF[lb].Kab = Kab;

      schF[la].PA[0] = schF[lb].PA[0] = PA[0];
      schF[la].PA[1] = schF[lb].PA[1] = PA[1];
      schF[la].PA[2] = schF[lb].PA[2] = PA[2];

      schF[la].PB[0] = schF[lb].PB[0] = PB[0];
      schF[la].PB[1] = schF[lb].PB[1] = PB[1];
      schF[la].PB[2] = schF[lb].PB[2] = PB[2];
    }
  }

}


void MatElements::IntegralEvaluation (int num, fstream& boysFile, fstream& boysFile1) {

  int n = (num * (num - 1) / 2);

  OrbitalData(num);
  int totAtm, gs, gsTS;
  totAtm = readData.lastBasis.size();
  atoms.resize(totAtm);

  bE.resize(n);
  bETS.resize(num);

  gs = ghF.checkGroupSize(n);
  gsTS = ghF.checkGroupSize(num);

//  cout << "voy con: " << n/gs << "   " << gs << std::endl;
//  cout << "voy con: " << num/gsTS << "   " << gsTS << std::endl;

  kineticEnergy.resize(readData.numberOfOccOrb);
  coulombicEnergy.resize(readData.numberOfOccOrb);
  oneElecCoulombEnergy.resize(readData.numberOfOccOrb);
  exchangeEnergy.resize(readData.numberOfOccOrb);

  double res;
  for (int i=0; i<totAtm; i++) {
    atoms[i].atmType = readData.lastBasis[i].atomType;
    atoms[i].coords[0] = readData.lastBasis[i].nucleusCoord[0];
    atoms[i].coords[1] = readData.lastBasis[i].nucleusCoord[1];
    atoms[i].coords[2] = readData.lastBasis[i].nucleusCoord[2];
  //  cout << atoms[i].atmType << "   " << atoms[i].coords[0] << "  " << atoms[i].coords[1] << "  " << atoms[i].coords[2] << std::endl;
  }

  LectureBoysValuesTest<double>(boysFile, boysFile1);

  queue Q{gpu_selector{}};

  {

  buffer<oneElectro<double>, 1> bufoE{oEi};
  buffer<oneElectro<double>, 1> bufoETS{oETS};
  buffer<bielectro<double>, 1> bufbETS{bETS};
  buffer<bielectro<double>, 1> bufbE{bE};

  buffer<schwarz<double>, 1> bufswF{schF.data(), schF.size()};
  buffer<atmInfo<double>, 1> bufAtms{atoms.data(), atoms.size()};

  buffer<double, 1> bufSCTNumA{btvNuma.data(), btvNuma.size()};
  buffer<double, 1> bufSCTNumB{btvNumb.data(), btvNumb.size()};
  buffer<double, 1> bufSCTNumC{btvNumc.data(), btvNumc.size()};
  buffer<double, 1> bufSCTNumD{btvNumd.data(), btvNumd.size()};

  buffer<double, 1> bufSCTA{btva.data(), btva.size()};
  buffer<double, 1> bufSCTB{btvb.data(), btvb.size()};
  buffer<double, 1> bufSCTC{btvc.data(), btvc.size()};
  buffer<double, 1> bufSCTD{btvd.data(), btvd.size()};

//////////////////////////////////////////////(i| j)////////////////////////////////////////////////////////////

  Q.submit([&](handler& h) {

    auto accoE = bufoE.get_access<access::mode::write>(h);
    auto accbE = bufbE.get_access<access::mode::write>(h);

    accessor accAtms{bufAtms, h, read_only};
    accessor accsw{bufswF, h, read_write};

    auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
    auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
    auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
    auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

    auto accSCA = bufSCTA.get_access<access::mode::read>(h);
    auto accSCB = bufSCTB.get_access<access::mode::read>(h);
    auto accSCC = bufSCTC.get_access<access::mode::read>(h);
    auto accSCD = bufSCTD.get_access<access::mode::read>(h);

    range num_groups = n/gs;
    range group_size = gs;

    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {

      storm<double> bF;
      int charge, angT, angs[3], ang1[3], ang2[3];
      double PA[3], PB[3]; 
      double Px, Py, Pz, factor, expo2, expoTot, Kab, traslapeTot, nucElecTot, kineticTot, Xpc, Ypc, Zpc, Rpc, nucleousElectron, integralNucEle, coefEx, coefEy, coefEz, kineticA, kineticB, kineticC, overlap[3], kinetic[3], coefs;
      double vala, valb, valc, vald, xi, res, boysInt, bx; 
      double boysDx[10], boysDy[10], boysDz[10], boysT[15]; 
      int ind, ind_j, local_i, EfactorsX[15], EfactorsY[15], EfactorsZ[15];
      int nx, ny, nz, nt, ni, i, j; 
      double val;
      int boysNx[10], boysNy[10], boysNz[10];
      auto ind_i = group_i*group_size + item.get_local_id();

      i = accoE[ind_i].ind_a; 
      j = accoE[ind_i].ind_b; 

      accbE[ind_i].ind_a = i;
      accbE[ind_i].ind_b = j;
      accbE[ind_i].ind_c = i;
      accbE[ind_i].ind_d = j;

      local_i = i * num + j;

      expoTot = accsw[local_i].expoTot;
      Kab = accsw[local_i].Kab;
      angT = accsw[local_i].angT;
      expo2 = accsw[local_i].expo2;

      Px = accsw[local_i].P[0];
      Py = accsw[local_i].P[1];
      Pz = accsw[local_i].P[2];

      PA[0] = accsw[local_i].PA[0];
      PA[1] = accsw[local_i].PA[1];
      PA[2] = accsw[local_i].PA[2];

      PB[0] = accsw[local_i].PB[0];
      PB[1] = accsw[local_i].PB[1];
      PB[2] = accsw[local_i].PB[2];

      angs[0] = accsw[local_i].angs[0];
      angs[1] = accsw[local_i].angs[1];
      angs[2] = accsw[local_i].angs[2];

      ang1[0] = accsw[local_i].ang1[0];
      ang1[1] = accsw[local_i].ang1[1];
      ang1[2] = accsw[local_i].ang1[2];

      ang2[0] = accsw[local_i].ang2[0];
      ang2[1] = accsw[local_i].ang2[1];
      ang2[2] = accsw[local_i].ang2[2];

      if (Kab < -40.0) {
        traslapeTot = 0.0;
        nucElecTot = 0.0;
        kineticTot = 0.0;
      } else {
        Kab = sycl::exp(Kab);
        nucElecTot = 0.0;
        nucleousElectron = 0.0;
        if (angT == 0) {
          for (int v=0; v<totAtm; v++) {
            charge = accAtms[v].atmType;
            Xpc = Px-accAtms[v].coords[0];
            Ypc = Py-accAtms[v].coords[1];
            Zpc = Pz-accAtms[v].coords[2];
            Rpc = Xpc*Xpc + Ypc*Ypc + Zpc*Zpc;

            factor = expoTot * Rpc;
            if (factor < 1E-12)
              factor = 0.0;

            nucleousElectron += charge * bF.boysZeroDegree<double>(factor);
          }
         nucElecTot = (-2.0*M_PI/expoTot)*Kab * nucleousElectron;
       } else {
         nucElecTot = 0.0;
        for (int k=0; k<totAtm; k++) {
          charge = accAtms[k].atmType;
          Xpc = Px-accAtms[k].coords[0];
          Ypc = Py-accAtms[k].coords[1];
          Zpc = Pz-accAtms[k].coords[2];
          Rpc = Xpc*Xpc + Ypc*Ypc + Zpc*Zpc;

          factor = expoTot * Rpc;

          integralNucEle = 0.0;
          if (factor < 1E-12)
            factor = 0;
          ind_j = factor * 10000;
          xi = factor - ind_j * 0.0001;

          boysT[0] = bF.boysZeroDegree<double>(factor);

          if (factor == 0 || factor >= 30.0) {
            for (int u=1; u<=angT; u++) 
              boysT[u] = bF.ObaraIntegral<double>(u, factor);
          } else { 
            nt = 3;
            if (angT <= nt)
              nt = angT;
            for (int i=1; i<=nt; i++) {
              ind = ind_j + (i-1)*300001;
              vala = accSCNumA[ind];
              valb = accSCNumB[ind];
              valc = accSCNumC[ind];
              vald = accSCNumD[ind];
              boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
              }
            if (angT > 3) {
              nt = angT / 5;
              if (nt > 0) {
                val = angT / 5.0;
                val -= nt;
                if (val == 0.0)
                nt--;
              }
              for (int u=0; u<=nt; u++) {
                ind = ind_j + u*300001;
                vala = accSCA[ind];
                valb = accSCB[ind];
                valc = accSCC[ind];
                vald = accSCD[ind];

                bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                if (u == 0)
                  bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                else
                  bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                }
              }
            }

          for (int t=0; t<=angs[0]; t++) 
            for (int u=0; u<=angs[1]; u++)
              for (int v=0; v<=angs[2]; v++) {
                bF.DerivativeBoys(boysDx, boysNx, expoTot, Xpc, t, nx);
                bF.DerivativeBoys(boysDy, boysNy, expoTot, Ypc, u, ny);
                bF.DerivativeBoys(boysDz, boysNz, expoTot, Zpc, v, nz);
                coefEx = bF.McMurchieCoefLoop(EfactorsX, ang1[0], ang2[0], t, expoTot, PA[0], PB[0]);
                coefEy = bF.McMurchieCoefLoop(EfactorsY, ang1[1], ang2[1], u, expoTot, PA[1], PB[1]);
                coefEz = bF.McMurchieCoefLoop(EfactorsZ, ang1[2], ang2[2], v, expoTot, PA[2], PB[2]);
                coefs = coefEx * coefEy * coefEz;

                integralNucEle = 0.0;
                for (int mx=0; mx<nx; mx++)
                  for (int my=0; my<ny; my++)
                    for (int mz=0; mz<nz; mz++) {
                      res = boysDx[mx]*boysDy[my]*boysDz[mz];
                      ni = boysNx[mx]+boysNy[my]+boysNz[mz];
                      integralNucEle += boysT[ni]*res;
                    }
                
                nucElecTot += integralNucEle * coefs * charge;
              }
          }
          nucElecTot *= (-2.0*M_PI/expoTot)*Kab;
        }
      traslapeTot = Kab;
      for (int i=0; i<3; i++) {
        overlap[i] = bF.overlapIntegral(ang1[i], ang2[i], PA[i], PB[i], expoTot);
        traslapeTot *= overlap[i];
        kineticA = -2.0*expo2*expo2*bF.overlapIntegral(ang1[i], ang2[i]+2, PA[i], PB[i], expoTot);
        kineticB = expo2*(2.0*ang2[i] + 1)*overlap[i];
        kineticC = -0.5*ang2[i]*(ang2[i]-1)*bF.overlapIntegral(ang1[i], ang2[i]-2, PA[i], PB[i], expoTot);
        kinetic[i] = kineticA + kineticB + kineticC;
       }

       kineticTot = (kinetic[0]*overlap[1]*overlap[2] + overlap[0]*kinetic[1]*overlap[2] + overlap[0]*overlap[1]*kinetic[2]);
      kineticTot *= Kab;
      }

    accoE[ind_i].nE = nucElecTot; 
    accoE[ind_i].kE = kineticTot;
    accoE[ind_i].oE = traslapeTot;
      
      });
    });
  });


//////////////////////////////////////////////(i | i)////////////////////////////////////////////////////////////

  Q.submit([&](handler& h) {

    auto accoE = bufoETS.get_access<access::mode::write>(h);
    auto accbE = bufbETS.get_access<access::mode::write>(h);

    accessor accAtms{bufAtms, h, read_only};
    accessor accsw{bufswF, h, read_write};

    auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
    auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
    auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
    auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

    auto accSCA = bufSCTA.get_access<access::mode::read>(h);
    auto accSCB = bufSCTB.get_access<access::mode::read>(h);
    auto accSCC = bufSCTC.get_access<access::mode::read>(h);
    auto accSCD = bufSCTD.get_access<access::mode::read>(h);

    range num_groups = num/gsTS;
    range group_size = gsTS;

    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id(0);
      grp.parallel_for_work_item([&](h_item<1> item) {

      storm<double> bF;
      int charge, angT, angs[3], ang1[3], ang2[3];
      double Px, Py, Pz, factor, expo2, expoTot, Kab, traslapeTot, nucElecTot, kineticTot, Xpc, Ypc, Zpc, Rpc, nucleousElectron, integralNucEle, coefEx, coefEy, coefEz, kineticA, kineticB, kineticC, overlap[3], kinetic[3], coefs;
      double vala, valb, valc, vald, xi, res, boysInt, bx;
      double boysDx[7], boysDy[7], boysDz[7], boysT[50];
      int ind, ind_j, local_i, i, j, EfactorsX[15], EfactorsY[15], EfactorsZ[15];
      int nx, ny, nz, nt, ni;
      double val;
      int boysNx[7], boysNy[7], boysNz[7];
      auto ind_i = group_i*group_size + item.get_local_id();

      i = accoE[ind_i].ind_a;

      accbE[ind_i].ind_a = i;
      accbE[ind_i].ind_b = i;
      accbE[ind_i].ind_c = i;
      accbE[ind_i].ind_d = i;
  
      local_i = i*num + i;

      expoTot = accsw[local_i].expoTot;
      Kab = accsw[local_i].Kab;
      angT = accsw[local_i].angT;
      expo2 = accsw[local_i].expo2;

      Px = accsw[local_i].P[0];
      Py = accsw[local_i].P[1];
      Pz = accsw[local_i].P[2];

      angs[0] = accsw[local_i].angs[0];
      angs[1] = accsw[local_i].angs[1];
      angs[2] = accsw[local_i].angs[2];

      ang1[0] = accsw[local_i].ang1[0];
      ang1[1] = accsw[local_i].ang1[1];
      ang1[2] = accsw[local_i].ang1[2];

      ang2[0] = accsw[local_i].ang2[0];
      ang2[1] = accsw[local_i].ang2[1];
      ang2[2] = accsw[local_i].ang2[2];

      Kab = 1.0;
      nucElecTot = 0.0;
      nucleousElectron = 0.0;
      if (angT == 0) {
        for (int v=0; v<totAtm; v++) {
          charge = accAtms[v].atmType;
          Xpc = Px-accAtms[v].coords[0];
          Ypc = Py-accAtms[v].coords[1];
          Zpc = Pz-accAtms[v].coords[2];
          Rpc = Xpc*Xpc + Ypc*Ypc + Zpc*Zpc;

          factor = expoTot * Rpc;
          if (factor < 1E-12)
            factor = 0.0;

          nucleousElectron += charge * bF.boysZeroDegree<double>(factor);
        }
       nucElecTot = (-2.0*M_PI/expoTot) * nucleousElectron;
      } else {
        nucElecTot = 0.0;
        for (int k=0; k<totAtm; k++) {
          charge = accAtms[k].atmType;
          Xpc = Px-accAtms[k].coords[0];
          Ypc = Py-accAtms[k].coords[1];
          Zpc = Pz-accAtms[k].coords[2];
          Rpc = Xpc*Xpc + Ypc*Ypc + Zpc*Zpc;

          factor = expoTot * Rpc;

          integralNucEle = 0.0;
          if (factor < 1E-12)
            factor = 0;
          ind_j = factor * 10000;
          xi = factor - ind_j * 0.0001;

          boysT[0] = bF.boysZeroDegree<double>(factor);

          if (factor == 0 || factor >= 30.0) {
            for (int u=1; u<=angT; u++)
              boysT[u] = bF.ObaraIntegral<double>(u, factor);
          } else {
            nt = 3;
            if (angT <= nt)
              nt = angT;
            for (int i=1; i<=nt; i++) {
              ind = ind_j + (i-1)*300001;
              vala = accSCNumA[ind];
              valb = accSCNumB[ind];
              valc = accSCNumC[ind];
              vald = accSCNumD[ind];
              boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
              }
            if (angT > 3) {
              nt = angT / 5;
              if (nt > 0) {
                val = angT / 5.0;
                val -= nt;
                if (val == 0.0)
                nt--;
              }
              for (int u=0; u<=nt; u++) {
                ind = ind_j + u*300001;
                vala = accSCA[ind];
                valb = accSCB[ind];
                valc = accSCC[ind];
                vald = accSCD[ind];

                bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                if (u == 0)
                  bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                else
                  bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                }
              }
            }

          for (int t=0; t<=angs[0]; t++)
            for (int u=0; u<=angs[1]; u++)
              for (int v=0; v<=angs[2]; v++) {
                bF.DerivativeBoys(boysDx, boysNx, expoTot, Xpc, t, nx);
                bF.DerivativeBoys(boysDy, boysNy, expoTot, Ypc, u, ny);
                bF.DerivativeBoys(boysDz, boysNz, expoTot, Zpc, v, nz);
                coefEx = bF.McMurchieCoefLoop(EfactorsX, ang1[0], ang2[0], t, expoTot, 0.0, 0.0);
                coefEy = bF.McMurchieCoefLoop(EfactorsY, ang1[1], ang2[1], u, expoTot, 0.0, 0.0);
                coefEz = bF.McMurchieCoefLoop(EfactorsZ, ang1[2], ang2[2], v, expoTot, 0.0, 0.0);
                coefs = coefEx * coefEy * coefEz;

                integralNucEle = 0.0;
                for (int mx=0; mx<nx; mx++)
                  for (int my=0; my<ny; my++)
                    for (int mz=0; mz<nz; mz++) {
                      res = boysDx[mx]*boysDy[my]*boysDz[mz];
                      ni = boysNx[mx]+boysNy[my]+boysNz[mz];
                      integralNucEle += boysT[ni]*res;
                    }

                nucElecTot += integralNucEle * coefs*charge;
              }
          }
          nucElecTot *= (-2.0*M_PI/expoTot);
        }
      traslapeTot = Kab;
      kineticTot = Kab;
      for (int i=0; i<3; i++) {
        overlap[i] = bF.overlapIntegral(ang1[i], ang2[i], 0.0, 0.0, expoTot);
        traslapeTot *= overlap[i];
        kineticA = -2.0*expo2*expo2*bF.otherOverlap(ang1[i], ang2[i]+2, 0.0, 0.0, expoTot, EfactorsX);
        kineticB = expo2*(2.0*ang2[i] + 1)*overlap[i];
        kineticC = -0.5*ang2[i]*(ang2[i]-1)*bF.otherOverlap(ang1[i], ang2[i]-2, 0.0, 0.0, expoTot, EfactorsY);
        kinetic[i] = kineticA + kineticB + kineticC;
       }

       kineticTot *= (kinetic[0]*overlap[1]*overlap[2] + overlap[0]*kinetic[1]*overlap[2] + overlap[0]*overlap[1]*kinetic[2]);

    accoE[ind_i].nE = nucElecTot;
    accoE[ind_i].kE = kineticTot;
    accoE[ind_i].oE = traslapeTot;

      });
    });
  });

///////////////////////////STARTS THE BIELECTRONIC PART ////////////////////////////////
///////////////////////////////////////(i i | i i)////////////////////////////////////////////////////////////

    Q.submit([&](handler& h) {

    auto accbE = bufbETS.get_access<access::mode::write>(h);

    accessor accswF{bufswF, h, read_write};

    range num_groups = num/gsTS;
    range group_size = gsTS;

    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, totAngMom, ang12[3], ang1[3];
        double expo12, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, scr;
        double boysDx[7], boysDy[7], boysDz[7], boysT[50];
        int ind, ind_j, local_i, EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int nx, ny, nz, ni, totalCenter, i;
        double factorCou, etaCou;
        int boysNx[7], boysNy[7], boysNz[7], la, lb;
        auto ind_i = group_i*group_size + item.get_local_id();

        i = accbE[ind_i].ind_a;
        local_i = i*num + i;

        angT12 = accswF[local_i].angT;

        expo12 = accswF[local_i].expoTot; 

        ang12[0] = accswF[local_i].angs[0];
        ang12[1] = accswF[local_i].angs[1];
        ang12[2] = accswF[local_i].angs[2];

        ang1[0] = accswF[local_i].ang1[0];
        ang1[1] = accswF[local_i].ang1[1];
        ang1[2] = accswF[local_i].ang1[2];


        factorCou = bF.pi52/(expo12*expo12*sycl::sqrt(2.0*expo12));
        etaCou = 0.5*expo12;
        totAngMom = 2*angT12;

        if (angT12 == 0)
          bielecTot = factorCou;
        else 
          bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang12, ang1, ang1, ang1, ang1, totAngMom, etaCou, expo12, expo12); 
     
        accbE[ind_i].val = bielecTot;

        if (bielecTot == 0.0)
          scr = 0.0;
        else 
          scr = sycl::sqrt(bielecTot);

        accswF[local_i].scr = scr;

        });
      });
    });

/////////////////////////////////////////////////////////(i j | i j)///////////////////////////////


    Q.submit([&](handler& h) {

    auto accbE = bufbE.get_access<access::mode::write>(h);

    accessor accsw{bufswF, h, read_write};

    range num_groups = n/gs;
    range group_size = gs;

    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
    grp.parallel_for_work_item([&](h_item<1> item) {
      storm<double> bF;
      int angT12, totAngMom, ang12[3], ang1[3], ang2[3];
      double expo12, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, totalCou, PQ[3], PA[3], PB[3];
      double boysDx[7], boysDy[7], boysDz[7], boysT[50], scr;
      int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
      int center1, center2, la, lb, i, j, local_i, local_j;
      double factorCou, etaCou;
      int boysNx[7], boysNy[7], boysNz[7];
      auto ind_i = group_i*group_size + item.get_local_id(); 

      i = accbE[ind_i].ind_a;
      j = accbE[ind_i].ind_b;

      local_i = i*num + j;
      local_j = j*num + i;
      
      angT12 = accsw[local_i].angT;

      expo12 = accsw[local_i].expoTot;

      ang12[0] = accsw[local_i].angs[0];
      ang12[1] = accsw[local_i].angs[1];
      ang12[2] = accsw[local_i].angs[2];

      ang1[0] = accsw[local_i].ang1[0];
      ang1[1] = accsw[local_i].ang1[1];
      ang1[2] = accsw[local_i].ang1[2];

      ang2[0] = accsw[local_i].ang2[0];
      ang2[1] = accsw[local_i].ang2[1];
      ang2[2] = accsw[local_i].ang2[2];

      factorCou = bF.pi52/(expo12*expo12*sycl::sqrt(2.0*expo12));
      etaCou = 0.5*expo12;
      totAngMom = angT12 + angT12;

      center1 = accsw[local_i].atmCenter1;
      center2 = accsw[local_i].atmCenter2;

      if (center1 == center2) {
        if (angT12 == 0)
          bielecTot = factorCou;
        else
          bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang12, ang1, ang2, ang1, ang2, totAngMom, etaCou, expo12, expo12);

      } else {
          Kab = accsw[local_i].Kab;
          if (Kab < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            totalCou = factorCou*Kab*Kab;

            PA[0] = accsw[local_i].PA[0];
            PA[1] = accsw[local_i].PA[1];
            PA[2] = accsw[local_i].PA[2];

            PB[0] = accsw[local_i].PB[0];
            PB[1] = accsw[local_i].PB[1];
            PB[2] = accsw[local_i].PB[2];

            PQ[0] = 0.0; 
            PQ[1] = 0.0; 
            PQ[2] = 0.0; 

            if (angT12 == 0) {
              bielecTot = 1.0;
            } else {
              boysT[0] = 1.0;
              for (int i=1; i<=totAngMom; i++)
                boysT[i] = bF.ObaraIntegral<double>(i, 0.0);

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang12, ang1, ang2, ang1, ang2, PQ, PA, PB, PA, PB, etaCou, expo12, expo12);
            }

            bielecTot *= totalCou;
          }
      }

      accbE[ind_i].val = bielecTot;

      if (bielecTot == 0.0)
        scr = 0.0;
      else
        scr = sycl::sqrt(bielecTot);

      accsw[local_i].scr = scr;
      accsw[local_j].scr = scr;

        });
      });
    });

  }Q.wait();

////////////////////////////////////////////////////VALUES PASSING///////////////////////////////////////////////////////////////////////////////////////////////////

  int n1, n2, n3, n4, n5, n6, n7, l1;
  int gs1, gs2, gs3, gs4, gs5, gs6, gs7;
////////(ii|ij)
  n1 = num*(num-1);
  bE1.resize(n1);
  gs1 = ghF.checkGroupSize(n1);

//  cout << "voy con: " << n1/gs1 << "   " << gs1 << std::endl;
  n1 = 0;
  for (int i=0; i<aO.size(); i++) 
    for (int j=0; j<aO.size(); j++) {
      if (j != i) {
        bE1[n1].ind_a = i;
        bE1[n1].ind_b = i;
        bE1[n1].ind_c = i;
        bE1[n1].ind_d = j;
        n1++;
      }
   }


/////////(ii|jj)
  bE2.resize(n);
  gs2 = gs;
//  cout << "voy con: " << n/gs2 << "   " << gs2 << std::endl;

  n2 = 0;
  int j1 = 1;
  for (int i=0; i<aO.size(); i++) {
    for (int j=j1; j<aO.size(); j++) {
      bE2[n2].ind_a = i;
      bE2[n2].ind_b = i;
      bE2[n2].ind_c = j;
      bE2[n2].ind_d = j;
      n2++;
    }
    j1++;
  }

//////(ii|jk)
  n3 = num*(num-1)*(num-2)/2;
  bE3.resize(n3);
  gs3 = ghF.checkGroupSize(n3);
  int k1;

  n3 = 0;
  for (int i=0; i<aO.size(); i++) {
    k1 = 1;
    for (int j=0; j<aO.size(); j++) {
      if (j == i)
        j++;
      for (int k=k1; k<aO.size(); k++) {
        if (k == j)
          k++;
        if (k == i)
          k++;
        if (k == j)
          k++;
        if (k < aO.size()) {
          bE3[n3].ind_a = i;
          bE3[n3].ind_b = i;
          bE3[n3].ind_c = j;
          bE3[n3].ind_d = k;
          n3++;
        }
      }
      k1++;
    }
  }


//////(0i|0j)
  n4 = (num-1)*(num-2)/2;
  bE4.resize(n4);
  gs4 = ghF.checkGroupSize(n4);
//  cout << "voy con: " << n4/gs4 << "   " << gs4 << std::endl;
  
  n4 = 0;
  for (int i=1; i<aO.size(); i++) 
    for (int j=i+1; j<aO.size(); j++) {
      bE4[n4].ind_a = 0;
      bE4[n4].ind_b = i;
      bE4[n4].ind_c = 0;
      bE4[n4].ind_d = j;
      n4++;
    }

//////(ij|ik)
  n5 = (num-1)*(num-1)*(num-2)/2;
  bE5.resize(n5);
  gs5 = ghF.checkGroupSize(n5);
//  cout << "voy con: " << n5/gs5 << "   " << gs5 << std::endl;

  n5 = 0;
  for (int i=1; i<aO.size(); i++) {
    for (int j=0; j<aO.size(); j++) {
      if (j == i)
        j++;
      k1 = j+1;
      for (int k=k1; k<aO.size(); k++) {
        if (k == j)
          k++;
        if (k == i)
          k++;
        if (k == j)
          k++;

        if (k < aO.size()) {
          bE5[n5].ind_a = i;
          bE5[n5].ind_b = j;
          bE5[n5].ind_c = i;
          bE5[n5].ind_d = k;
          n5++;
        }
      }
      k1++;
      }
    }
  
///////(0i|jk)
  n6 = (num-1)*(num-2)*(num-3)/2;
  bE6.resize(n6);
  gs6 = ghF.checkGroupSize(n6);
//  cout << "voy con: " << n6/gs6 << "   " << gs6 << std::endl;

  n6 = 0;
  for (int j=1; j<aO.size(); j++) {
    k1 = 1;
    for (int k=k1; k<aO.size(); k++) {
      if (k == j)
        k++;
      l1 = k + 1;
      for (int l=l1; l<aO.size(); l++) {
        if (l == j)
          l++;
        if (l < aO.size()) {
          bE6[n6].ind_a = 0;
          bE6[n6].ind_b = j;
          bE6[n6].ind_c = k;
          bE6[n6].ind_d = l;
          n6++;
        }
      }
      k1++;
    }
  }

//////////(ij|kl)
  n7 = (num-1)*(num-2)*(num-3)*(num-4)/8;
  bE7.resize(n7);
  gs7 = ghF.checkGroupSize(n7);
//  cout << "voy con: " << n7/gs7 << "   " << gs7 << std::endl;

  n7 = 0;
  for (int i=1; i<aO.size(); i++) {
    j1 = i+1;
    k1 = i;
    for (int j=j1; j<aO.size(); j++) {
      for (int k=k1; k<aO.size(); k++) {
        if (k == i)
          k++;
        if (k == j)
          k++;

        l1 = k+1;
        for (int l=l1; l<aO.size(); l++) {
        if (l == j)
          l++;
        if (l < aO.size()) {
          bE7[n7].ind_a = i;
          bE7[n7].ind_b = j;
          bE7[n7].ind_c = k;
          bE7[n7].ind_d = l;
          n7++;
          }
        }
      }
      j1++;
    }
  }



////////////////////////////////////////////////////////(i i | i j) ////////////////////////////////////////////////////////////////////////////////////////////////

  queue W{gpu_selector{}};

  {

  buffer<schwarz<double>, 1> bufswF{schF.data(), schF.size()};

  buffer<double, 1> bufSCTNumA{btvNuma.data(), btvNuma.size()};
  buffer<double, 1> bufSCTNumB{btvNumb.data(), btvNumb.size()};
  buffer<double, 1> bufSCTNumC{btvNumc.data(), btvNumc.size()};
  buffer<double, 1> bufSCTNumD{btvNumd.data(), btvNumd.size()};

  buffer<double, 1> bufSCTA{btva.data(), btva.size()};
  buffer<double, 1> bufSCTB{btvb.data(), btvb.size()};
  buffer<double, 1> bufSCTC{btvc.data(), btvc.size()};
  buffer<double, 1> bufSCTD{btvd.data(), btvd.size()};

  buffer<bielectro<double>, 1> bufbE1{bE1};
  buffer<bielectro<double>, 1> bufbE2{bE2};
  buffer<bielectro<double>, 1> bufbE3{bE3};
  buffer<bielectro<double>, 1> bufbE4{bE4};
  buffer<bielectro<double>, 1> bufbE5{bE5};
  buffer<bielectro<double>, 1> bufbE6{bE6};
  buffer<bielectro<double>, 1> bufbE7{bE7};

   W.submit([&](handler& h) {

   auto accbE = bufbE1.get_access<access::mode::read_write>(h);

   accessor accsw{bufswF, h, read_only};
 
   auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
   auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
   auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
   auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

   auto accSCA = bufSCTA.get_access<access::mode::read>(h);
   auto accSCB = bufSCTB.get_access<access::mode::read>(h);
   auto accSCC = bufSCTC.get_access<access::mode::read>(h);
   auto accSCD = bufSCTD.get_access<access::mode::read>(h);

   range num_groups = n1/gs1;
   range group_size = gs1;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();
 
      ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
      ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;
     
      expo12 = accsw[ia].expoTot;
      expo34 = accsw[ja].expoTot;
  
      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];   

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne <= bF.tolerance)
        bielecTot = 0.0;
      else { 
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;
  
        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;
       
        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else 
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];
 
              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });


///////////////////////////////////////////////////////(i i | j j)////////////////////////////////////////////////////////////////////////////////////////////


   W.submit([&](handler& h) {

   auto accbE = bufbE2.get_access<access::mode::read_write>(h);

   accessor accsw{bufswF, h, read_only};

   auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
   auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
   auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
   auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

   auto accSCA = bufSCTA.get_access<access::mode::read>(h);
   auto accSCB = bufSCTB.get_access<access::mode::read>(h);
   auto accSCC = bufSCTC.get_access<access::mode::read>(h);
   auto accSCD = bufSCTD.get_access<access::mode::read>(h);

   range num_groups = n2/gs2;
   range group_size = gs2;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();

        ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
        ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;

        expo12 = accsw[ia].expoTot;
        expo34 = accsw[ja].expoTot;

      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne <= bF.tolerance)
        bielecTot = 0.0;
      else {
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;

        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;

        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];

              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });



/////////////////////////////////////////////////////////(i i| j k)////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  W.submit([&](handler& h) {

    auto accbE = bufbE3.get_access<access::mode::read_write>(h);

    accessor accsw{bufswF, h, read_only};

    auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
    auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
    auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
    auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

    auto accSCA = bufSCTA.get_access<access::mode::read>(h);
    auto accSCB = bufSCTB.get_access<access::mode::read>(h);
    auto accSCC = bufSCTC.get_access<access::mode::read>(h);
    auto accSCD = bufSCTD.get_access<access::mode::read>(h);

    range num_groups = n3/gs3;
    range group_size = gs3;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();

      ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
      ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;

      expo12 = accsw[ia].expoTot;
      expo34 = accsw[ja].expoTot;

      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne < bF.tolerance)
        bielecTot = 0.0;
      else {
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;

        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;

        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];

              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });


///////////////////////////////////////////////////////(0 i| 0 j)/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   W.submit([&](handler& h) {

   auto accbE = bufbE4.get_access<access::mode::read_write>(h);

   accessor accsw{bufswF, h, read_only};

   auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
   auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
   auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
   auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

   auto accSCA = bufSCTA.get_access<access::mode::read>(h);
   auto accSCB = bufSCTB.get_access<access::mode::read>(h);
   auto accSCC = bufSCTC.get_access<access::mode::read>(h);
   auto accSCD = bufSCTD.get_access<access::mode::read>(h);

   range num_groups = n4/gs4;
   range group_size = gs4;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();

      ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
      ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;

      expo12 = accsw[ia].expoTot;
      expo34 = accsw[ja].expoTot;

      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne < bF.tolerance)
        bielecTot = 0.0;
      else {
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;

        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;

        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];

              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });


///////////////////////////////////////////////////(i j | i k)////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   W.submit([&](handler& h) {

   auto accbE = bufbE5.get_access<access::mode::read_write>(h);

   accessor accsw{bufswF, h, read_only};

   auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
   auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
   auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
   auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

   auto accSCA = bufSCTA.get_access<access::mode::read>(h);
   auto accSCB = bufSCTB.get_access<access::mode::read>(h);
   auto accSCC = bufSCTC.get_access<access::mode::read>(h);
   auto accSCD = bufSCTD.get_access<access::mode::read>(h);

   range num_groups = n5/gs5;
   range group_size = gs5;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();

      ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
      ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;

      expo12 = accsw[ia].expoTot;
      expo34 = accsw[ja].expoTot;

      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne <= bF.tolerance)
        bielecTot = 0.0;
      else {
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;

        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;

        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];

              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });


//////////////////////////////////////////////////(0 j | k l)//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  W.submit([&](handler& h) {

    auto accbE = bufbE6.get_access<access::mode::read_write>(h);

    accessor accsw{bufswF, h, read_only};

    auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
   auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
   auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
   auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

   auto accSCA = bufSCTA.get_access<access::mode::read>(h);
   auto accSCB = bufSCTB.get_access<access::mode::read>(h);
   auto accSCC = bufSCTC.get_access<access::mode::read>(h);
   auto accSCD = bufSCTD.get_access<access::mode::read>(h);

   range num_groups = n6/gs6;
   range group_size = gs6;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();

      ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
      ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;

      expo12 = accsw[ia].expoTot;
      expo34 = accsw[ja].expoTot;

      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne <= bF.tolerance)
        bielecTot = 0.0;
      else {
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;

        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;

        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];

              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });

////////////////////////////////////////////////////////////((i j | k l)//////////////////////////////////////////////////////////////////////////////////

  W.submit([&](handler& h) {

    auto accbE = bufbE7.get_access<access::mode::read_write>(h);

    accessor accsw{bufswF, h, read_only};

    auto accSCNumA = bufSCTNumA.get_access<access::mode::read>(h);
    auto accSCNumB = bufSCTNumB.get_access<access::mode::read>(h);
    auto accSCNumC = bufSCTNumC.get_access<access::mode::read>(h);
    auto accSCNumD = bufSCTNumD.get_access<access::mode::read>(h);

    auto accSCA = bufSCTA.get_access<access::mode::read>(h);
    auto accSCB = bufSCTB.get_access<access::mode::read>(h);
    auto accSCC = bufSCTC.get_access<access::mode::read>(h);
    auto accSCD = bufSCTD.get_access<access::mode::read>(h);

    range num_groups = n7/gs7;
    range group_size = gs7;
    h.parallel_for_work_group(num_groups, group_size, [=](group<1> grp) {
      int group_i = grp.get_id();
      grp.parallel_for_work_item([&](h_item<1> item) {
        storm<double> bF;
        int angT12, angT34, totAngMom, ang12[3], ang34[3], ang1[3], ang2[3], ang3[3], ang4[3];
        double expo12, expo34, bielecTot, coef1Ex, coef1Ey, coef1Ez, coefs1, coef2Ex, coef2Ey, coef2Ez, coefs2, integralBie, res, Kab, Kcd, totalCou, PQ[3], PA[3], PB[3], QC[3], QD[3], P[3], Q[3], Rpq, schwarzIne;
        double boysDx[10], boysDy[10], boysDz[10], boysT[35];
        int EfactorsX[15], EfactorsY[15], EfactorsZ[15];
        int center1, center2, center3, center4, ia, ja, ind, ind_j, nt, nx, ny, nz, ni;
        double vala, valb, valc, vald, val, xi, bx;
        double factorCou, etaCou, factor;
        int boysNx[10], boysNy[10], boysNz[10];

        auto local_i = group_i*group_size + item.get_local_id();

      ia = num*accbE[local_i].ind_a + accbE[local_i].ind_b;
      ja = num*accbE[local_i].ind_c + accbE[local_i].ind_d;

      expo12 = accsw[ia].expoTot;
      expo34 = accsw[ja].expoTot;

      P[0] = accsw[ia].P[0];
      P[1] = accsw[ia].P[1];
      P[2] = accsw[ia].P[2];

      Q[0] = accsw[ja].P[0];
      Q[1] = accsw[ja].P[1];
      Q[2] = accsw[ja].P[2];

      PQ[0] = P[0] - Q[0];
      PQ[1] = P[1] - Q[1];
      PQ[2] = P[2] - Q[2];

      Rpq = 0.0;
      for (int i=0; i<3; i++)
        Rpq += PQ[i]*PQ[i];

      schwarzIne = sycl::abs(accsw[ia].scr*accsw[ja].scr/(sycl::sqrt(Rpq) - accsw[ia].ext-accsw[ja].ext));

      if (schwarzIne <= bF.tolerance)
        bielecTot = 0.0;
      else {
        angT12 = accsw[ia].angT;
        angT34 = accsw[ja].angT;
        totAngMom = angT12 + angT34;

        ang12[0] = accsw[ia].angs[0];
        ang12[1] = accsw[ia].angs[1];
        ang12[2] = accsw[ia].angs[2];

        ang1[0] = accsw[ia].ang1[0];
        ang1[1] = accsw[ia].ang1[1];
        ang1[2] = accsw[ia].ang1[2];

        ang2[0] = accsw[ia].ang2[0];
        ang2[1] = accsw[ia].ang2[1];
        ang2[2] = accsw[ia].ang2[2];

        ang34[0] = accsw[ja].angs[0];
        ang34[1] = accsw[ja].angs[1];
        ang34[2] = accsw[ja].angs[2];

        ang3[0] = accsw[ja].ang1[0];
        ang3[1] = accsw[ja].ang1[1];
        ang3[2] = accsw[ja].ang1[2];

        ang4[0] = accsw[ja].ang2[0];
        ang4[1] = accsw[ja].ang2[1];
        ang4[2] = accsw[ja].ang2[2];

        factorCou = bF.pi52/(expo12*expo34*sycl::sqrt(expo12 + expo34));
        etaCou = expo12*expo34/(expo12+expo34);

        center1 = accsw[ia].atmCenter1;
        center2 = accsw[ia].atmCenter2;
        center3 = accsw[ja].atmCenter1;
        center4 = accsw[ja].atmCenter2;

        if (center1 == center2 && center3 == center4 && center1 == center3) {
          if (totAngMom == 0)
            bielecTot = factorCou;
          else
            bielecTot = factorCou*bF.OneCenterNewIntegral<double>(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, totAngMom, etaCou, expo12, expo34);
        } else {
          Kab = accsw[ia].Kab;
          Kcd = accsw[ja].Kab;
          if (Kab < -40.0 || Kcd < -40.0)
            bielecTot = 0.0;
          else {
            Kab = sycl::exp(Kab);
            Kcd = sycl::exp(Kcd);
            totalCou = factorCou*Kab*Kcd;
            factor = Rpq * etaCou;
            if (factor < 1E-12)
              factor = 0.0;
            if (totAngMom == 0)
              bielecTot = totalCou * bF.boysZeroDegree<double>(factor);
            else {

              PA[0] = accsw[ia].PA[0];
              PA[1] = accsw[ia].PA[1];
              PA[2] = accsw[ia].PA[2];

              PB[0] = accsw[ia].PB[0];
              PB[1] = accsw[ia].PB[1];
              PB[2] = accsw[ia].PB[2];

              QC[0] = accsw[ja].PA[0];
              QC[1] = accsw[ja].PA[1];
              QC[2] = accsw[ja].PA[2];

              QD[0] = accsw[ja].PB[0];
              QD[1] = accsw[ja].PB[1];
              QD[2] = accsw[ja].PB[2];

              ind_j = factor * 10000;
              xi = factor - ind_j * 0.0001;

              boysT[0] = bF.boysZeroDegree<double>(factor);

              if (factor == 0 || factor >= 30.0) {
                for (int u=1; u<=totAngMom; u++)
                  boysT[u] = bF.ObaraIntegral<double>(u, factor);
              } else {
                nt = 3;
                if (totAngMom <= nt)
                  nt = totAngMom;
              for (int i=1; i<=nt; i++) {
                ind = ind_j + (i-1)*300001;
                vala = accSCNumA[ind];
                valb = accSCNumB[ind];
                valc = accSCNumC[ind];
                vald = accSCNumD[ind];
                boysT[i] = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                }
              if (totAngMom > 3) {
                nt = totAngMom / 5;
                if (nt > 0) {
                  val = totAngMom / 5.0;
                  val -= nt;
                  if (val == 0.0)
                  nt--;
                }
                for (int u=0; u<=nt; u++) {
                  ind = ind_j + u*300001;
                  vala = accSCA[ind];
                  valb = accSCB[ind];
                  valc = accSCC[ind];
                  vald = accSCD[ind];

                  bx = bF.SplinesBoys<double>(xi, vala, valb, valc, vald);
                  if (u == 0)
                    bF.boysDownward(factor, bx, (u+1)*5, 4, boysT);
                  else
                    bF.boysDownward(factor, bx, (u+1)*5, u*5+1, boysT);
                  }
                }
              }

              bielecTot = bF.NewIntegral(boysT, boysDx, boysDy, boysDz, boysNx, boysNy, boysNz, EfactorsX, EfactorsY, EfactorsZ, ang12, ang34, ang1, ang2, ang3, ang4, PQ, PA, PB, QC, QD, etaCou, expo12, expo34);

              bielecTot *= totalCou;
            }
          }
        }
    }
     accbE[local_i].val = bielecTot;

       });
     });
   });


///////////////////////////////////////////////FINISHED//////////////////////////////////////////////////////////////////////////////

 }W.wait();


}

//////////////////////////////////////////////////////////////////////////////////////////


void MatElements::detectDevice () {
    // Loop through available devices in this platform

  for (auto const& this_platform : platform::get_platforms() ) {
    std::cout << "Found platform: "
    << this_platform.get_info<info::platform::name>() << "\n";

    for (auto const& Q : this_platform.get_devices() ) {
      std::cout << " Device: "
        << Q.get_info<info::device::name>() << "\n";
    }
    std::cout << "\n";
  }

  for (auto const& this_platform : platform::get_platforms() ) {
    std::cout << "Found Platform:\n";
    do_query<info::platform::name>(this_platform, "info::platform::name");
    do_query<info::platform::vendor>(this_platform, "info::platform::vendor");
    do_query<info::platform::version>(this_platform, "info::platform::version");
    do_query<info::platform::profile>(this_platform, "info::platform::profile");
// Loop through the devices available in this plaform
    for (auto &dev : this_platform.get_devices() ) {
      std::cout << " Device: "
        << dev.get_info<info::device::name>() << "\n";
      std::cout << " is_host(): "
        << (dev.is_host() ? "Yes" : "No") << "\n";
      std::cout << " is_cpu(): "
        << (dev.is_cpu() ? "Yes" : "No") << "\n";
      std::cout << " is_gpu(): "
        << (dev.is_gpu() ? "Yes" : "No") << "\n";
      std::cout << " is_accelerator(): "
        << (dev.is_accelerator() ? "Yes" : "No") << "\n";
      do_query<info::device::vendor>(dev, "info::device::vendor");
      do_query<info::device::driver_version>(dev, "info::device::driver_version");
      do_query<info::device::max_work_item_dimensions>(dev, "info::device::max_work_item_dimensions");
      do_query<info::device::max_work_group_size>(dev, "info::device::max_work_group_size");
      do_query<info::device::global_mem_size>(dev, "info::device::global_mem_size");
      do_query<info::device::local_mem_size>(dev, "info::device::local_mem_size");
      do_query<info::device::mem_base_addr_align>(dev, "info::device::mem_base_addr_align");
      do_query<info::device::partition_max_sub_devices>(dev, "info::device::partition_max_sub_devices");
      std::cout << " Many more queries are available than shown here!\n";
    }
    std::cout << "\n";
  }

}

void MatElements::Data(fstream& basisFile, fstream& xyzFile, fstream& movecsFile, fstream& boysFile, fstream& boysFile1) {

  readData.coefficientsFinal(basisFile, xyzFile, movecsFile);
  readData.nucleusBasis();
  readData.densityMatrix();
//  readData.getNumOccOrb();
  fourIndexTransf();
  
}

void MatElements::matrixElements (fstream& boysFile, fstream& boysFile1) {

	int num;

  num = 0;
  for (int i=0; i<readData.lastBasis.size(); i++)
    num += readData.lastBasis[i].basisExponents.size();

  IntegralEvaluation(num, boysFile, boysFile1);

}

void MatElements::matrixBuilding() { 

	int num;
  int ind_i, ind_j;

  num = 0;
  for (int i=0; i<readData.lastBasis.size(); i++)
    num += readData.lastBasis[i].basisExponents.size();

  cE.resize(num*num*num*num);
  oE.resize(num*num);
  kE.resize(num*num);
  nE.resize(num*num);
  
  for (int i=0; i<oEi.size(); i++) {
    ind_i = oEi[i].ind_a*num + oEi[i].ind_b;
    ind_j = oEi[i].ind_b*num + oEi[i].ind_a;
    oE[ind_i] = oEi[i].oE;
    kE[ind_i] = oEi[i].kE;
    nE[ind_i] = oEi[i].nE;
    oE[ind_j] = oEi[i].oE;
    kE[ind_j] = oEi[i].kE;
    nE[ind_j] = oEi[i].nE;
  }

  for (int i=0; i<oETS.size(); i++) {
    ind_i = oETS[i].ind_a*num + oETS[i].ind_b;
    oE[ind_i] = oETS[i].oE;
    kE[ind_i] = oETS[i].kE;
    nE[ind_i] = oETS[i].nE;
  }

  for (int i=0; i<bETS.size(); i++) {
    ind_i = bETS[i].ind_a;
    ind_j = num*num*num*ind_i + num*num*ind_i + num*ind_i + ind_i;
    cE[ind_j] = bETS[i].val;   
  }

  passingElements(num, bE, cE);
  passingElements(num, bE1, cE);
  passingElements(num, bE2, cE);
  passingElements(num, bE3, cE);
  passingElements(num, bE4, cE);
  passingElements(num, bE5, cE);
  passingElements(num, bE6, cE);
  passingElements(num, bE7, cE);

  int n, h;

  h = 0;
  n = readData.numberOfPrimFunc;
  // Aqui se arma la matriz Coulombica
  totalCoulombicEnergy = 0.0;
  double totalExchangeEnergy = 0.0;
  for (int i=0; i<readData.numberOfOccOrb; i++)  
    for (int j=0; j<readData.numberOfOccOrb; j++)  
      for (int a=0; a<readData.numberOfPrimFunc; a++)
        for (int b=0; b<readData.numberOfPrimFunc; b++) 
          for (int c=0; c<readData.numberOfPrimFunc; c++)
            for (int d=0; d<readData.numberOfPrimFunc; d++) {
              totalCoulombicEnergy += readData.dE[i*n*n+n*a+b]  * readData.dE[j*n*n+d*n+c] *cE[n*n*n*a + n*n*b + n*c + d];
              totalExchangeEnergy += readData.dE[i*n*n+n*a+b]  * readData.dE[j*n*n+d*n+c] *(-0.5*cE[n*n*n*a + n*n*c + n*b + d]);
              h++;
            }

        
  totalKineticEnergy = 0.0;
  totalNucleousElecEnergy = 0.0;
  for (int i=0; i<readData.numberOfOccOrb; i++)  
    for (int a=0; a<readData.numberOfPrimFunc; a++)
      for (int b=0; b<readData.numberOfPrimFunc; b++) {
        totalKineticEnergy += readData.dE[i*n*n+n*a+b] * kE[n*a+b];
        totalNucleousElecEnergy += readData.dE[i*n*n+n*a+b] * nE[n*a+b];
      }

//-------------------------------------TO PRINT OVERLAP MATRIX------------------

  overlapMatrix(num);

//------------------------------------------------------------------------------

  cout << "Bielectronic contribution:   " << 0.5*(totalCoulombicEnergy + totalExchangeEnergy) << std::endl;
  cout << "Monoelectronic contribution:   " << totalNucleousElecEnergy + totalKineticEnergy << std::endl;
  cout << "Energy Total   " << totalNucleousElecEnergy + totalKineticEnergy + 0.5*(totalCoulombicEnergy+totalExchangeEnergy) << std::endl;

}

void MatElements::overlapMatrix(int num) {

  int cont, cont1, tot, ni;

  vector<vector<double>> overlapElements;

  tot = readData.numberOfPrimFunc * readData.numberOfPrimFunc;

  cont = 0;
  for (int j=0; j<readData.lastBasis.size(); j++)
    for (int k=0; k<readData.lastBasis[j].shellSize.size(); k++)
      cont++;

  overlapElements.resize(cont);
  for (int i=0; i<cont; i++) {
    overlapElements[i].resize(cont);
    for (int j=0; j<cont; j++)
      overlapElements[i][j] = 0.0;
  }

  ni = cont;

  int a, b, ia, ib;

  cont1 = 0;
  b = 0;
  for (int j=0; j<readData.lastBasis.size(); j++) {
    ib = 0;
    for (int k=0; k<readData.lastBasis[j].shellSize.size(); k++) {
      for (int l=0; l<readData.lastBasis[j].shellSize[k]; l++) {
        cont = 0;
        a = 0;
        for (int h=0; h<readData.lastBasis.size(); h++) {
          ia = 0;
          for (int m=0; m<readData.lastBasis[h].shellSize.size(); m++) {
            for (int n=0; n<readData.lastBasis[h].shellSize[m]; n++) {
              overlapElements[b][a] += readData.lastBasis[j].basisCoefficients[ib] * readData.lastBasis[h].basisCoefficients[ia] * oE[num*cont1+cont];
              cont++;
              ia++;
            }
           a++;
          }
        }
        cont1++;
      ib++;
      }
      b++;
      }
    }

  for (int i=0; i<ni; i++) {
    for (int j=0; j<ni; j++) {
      if (i == j)
        cout << i << "  "  << j << "   " << overlapElements[i][j] << "\n";
    }
    cout << "---------------------------------------------------------------\n";
  }
}

void MatElements::passingElements(int n, vector<bielectro<double>> bE, vector<double>& cE) {
  
  int i, j, k, l;
  int n1, n2, n3, n4, n5, n6, n7, n8;
  for (int ia=0; ia<bE.size(); ia++) {
    i = bE[ia].ind_a;
    j = bE[ia].ind_b;
    k = bE[ia].ind_c;
    l = bE[ia].ind_d;
    n1 = n*n*n*i + n*n*j + n*k + l; // (ij|kl)
    n2 = n*n*n*i + n*n*j + n*l + k; // (ij|lk)
    n3 = n*n*n*j + n*n*i + n*k + l; // (ji|kl)
    n4 = n*n*n*j + n*n*i + n*l + k; // (ji|lk)
    n5 = n*n*n*l + n*n*k + n*i + j; // (lk|ij)
    n6 = n*n*n*l + n*n*k + n*j + i; // (lk|ji)
    n7 = n*n*n*k + n*n*l + n*i + j; // (kl|ij)
    n8 = n*n*n*k + n*n*l + n*j + i; // (kl|ji)
    cE[n1] = cE[n2] = cE[n3] = cE[n4] = cE[n5] = cE[n6] = cE[n7] = cE[n8] = bE[ia].val;
  }

}


void MatElements::fourIndexTransf() {

  int a;
  aO.resize(readData.numberOfPrimFunc);
  a = 0;
  for (int k=0; k<readData.lastBasis.size(); k++)
    for (int n=0; n<readData.lastBasis[k].basisExponents.size(); n++) {
      aO[a].atomCenter = k;
      aO[a].expo = readData.lastBasis[k].basisExponents[n];
      aO[a].funcType = readData.lastBasis[k].functionType[n];
      aO[a].angMom = readData.lastBasis[k].angularMoment[n];
      basis.angularMoment(aO[a].ang, aO[a].angMom, aO[a].funcType);

      aO[a].coord[0] = readData.lastBasis[k].nucleusCoord[0];
      aO[a].coord[1] = readData.lastBasis[k].nucleusCoord[1];
      aO[a].coord[2] = readData.lastBasis[k].nucleusCoord[2];
      a++;
    }


}

/*
void MatElements::MatrixBuilding() {

  size_t n = readData.numberOfPrimFunc;
  size_t occ = readData.numberOfOccOrb;

  repulsionElecEnergy.resize(occ*occ);
  totalNucleousElecEnergy = 0.0;
  totalKineticEnergy = 0.0;
  totalCoulombicEnergy = 0.0;

  {
    queue mM{gpu_selector{}};
    
    auto totalNucleousEnergy_dev = sycl::malloc_device<double>(1, mM);
    auto totalKineticEnergy_dev = sycl::malloc_device<double>(1, mM);
    auto totalCoulombicEnergy_dev = sycl::malloc_device<double>(1, mM);

    auto e1 = mM.memcpy(totalNucleousEnergy_dev, &totalNucleousElecEnergy, sizeof(double));
    auto e2 = mM.memcpy(totalKineticEnergy_dev, &totalKineticEnergy, sizeof(double));
    auto e3 = mM.memcpy(totalCoulombicEnergy_dev, &totalCoulombicEnergy, sizeof(double));

    buffer dM{readData.dE};
    buffer kM{kE};
    buffer nM{nE};
    buffer nEE{nucleousElecEnergy};
    buffer kineEne{kineticEnergy};
    buffer cM{bEt};
    buffer reEM{repulsionElecEnergy};

    mM.submit([&](handler &h) {
      accessor accdM(dM, h, read_only);
      accessor acckM(kM, h, read_only);
      accessor accnM(nM, h, read_only);
      accessor acccM(cM, h, read_only);
      accessor accnEE(nEE, h, write_only);
      accessor acckineEne(kineEne, h, write_only);
      accessor accrEE(reEM, h, write_only);
      h.depends_on(e1);
      h.depends_on(e2);
      h.depends_on(e3);
      h.single_task([=] () {
        for (int u=0; u<occ; u++) {
          acckineEne[u] = 0.0;
          accnEE[u] = 0.0;
          for (int i=0; i<n; i++)
            for (int j=0; j<n; j++) {
              acckineEne[u] += accdM[u*n*n+n*i+j] * acckM[n*i+j];
              accnEE[u] += accdM[u*n*n+n*i+j] * accnM[n*i+j];
            }
            (*totalNucleousEnergy_dev) += accnEE[u];
            (*totalKineticEnergy_dev) += acckineEne[u];
          for (int j=0; j<occ; j++) {
            for (int a=0; a<n; a++)
              for (int b=0; b<n; b++)
                for (int c=0; c<n; c++)
                  for (int d=0; d<n; d++) {
                    accrEE[u*occ+j] += accdM[u*n*n+n*a+b]  * accdM[j*n*n+c*n+d] * (acccM[n*n*n*a + n*n*b + n*c + d] - 0.5*acccM[n*n*n*a + n*n*c + n*b + d]);
                  }
            (*totalCoulombicEnergy_dev) += accrEE[u*occ+j];
            }
          }
      });
    });

  mM.memcpy(&totalNucleousElecEnergy, totalNucleousEnergy_dev, sizeof(double)).wait();
  mM.memcpy(&totalKineticEnergy, totalKineticEnergy_dev, sizeof(double)).wait();
  mM.memcpy(&totalCoulombicEnergy, totalCoulombicEnergy_dev, sizeof(double)).wait();

  sycl::free(totalNucleousEnergy_dev, mM);
  sycl::free(totalKineticEnergy_dev, mM);
  sycl::free(totalCoulombicEnergy_dev, mM);

  }

  cout << "Monoelectronic contribution:  " << totalNucleousElecEnergy + totalKineticEnergy << "\n";
  cout << "Bielectronic contribution:   " << 0.5*totalCoulombicEnergy << "\n";
  cout << "Energy Total   " << totalNucleousElecEnergy + totalKineticEnergy + 0.5*totalCoulombicEnergy<< "\n";

}
*/

