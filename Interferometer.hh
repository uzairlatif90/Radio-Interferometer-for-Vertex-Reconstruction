#ifndef INTER_HEAD
#define INTER_HEAD

#include "IceRayTracing.cc"
#include "gsl/gsl_multimin.h"
#include "AntennaConfig.hh"
#include "TRandom3.h"

using namespace std;

namespace Interferometer{

  void GenerateChHitTime(double DummyTx[3],double ChHitTime[2][TotalAntennasRx]);

  void ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]);

  void FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx]);

  void AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]);

  double MinimizeXYZ(const gsl_vector *v, void *params);
  
  void DoTheMinimization(double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx]);
 
}
#endif
