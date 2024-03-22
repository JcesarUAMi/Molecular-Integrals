#include <iostream>
#include <CL/sycl.hpp>
#include <fstream>
#include <cmath>
#include "BoysFunction.h"
using namespace std;
using namespace sycl;


template <typename U> class storm {
  public:
    template<typename T> T boysZeroDegree(T x) {
      T result, sqX;
      double sqPi = 0.5*sqrt(M_PI);

      if (x == 0.0)
        result = 1.0;
      else {
        sqX = sqrt(x);
        if (sqX >= 6.0)
          result = sqPi/sqX;
        else
          result = erf(sqX)*sqPi/sqX;
      }
      return result;
    }
};

template <typename U> class McMurchie {
  public:
    template <typename T> void EvaluationBoys (const vector<T>& a, T value) {
      buffer<T> buffA{a};
      
      queue Q;

      Q.submit([&](handler& h) {
        accessor matA{buffA, h};
        h.single_task([=]() {
          double result;
          storm<double> bF;  
          result = bF.boysZeroDegree(value);
          });
      });

    Q.wait();
    
    }


};

int main () {

  fstream boys;
  fstream boys1;

  double resultHost, resultHost1;
  double factor = 4.1346549;
 
  vector<double> ijole; 
  ijole.resize(5000000);

  double b[3], c[3];

  vector<double> a;
  a.resize(3);

  for (int i=0; i<3; i++) {
    a[i] = i*i;
    b[i] = 3*i;
    c[i] = 0;
  }


//  Boys<double> boysF;
//  resultHost = 0.0;
//  boysF.LectureBoysValuesTest<double>(boys, boys1);


//  resultHost = boysF.boysZeroDegree<double>(factor);
//  resultHost1 = boysF.BoysFunction<double>(2, factor);
   

  McMurchie<double> MC;
  MC.EvaluationBoys(a, factor);


  return 0;
}
