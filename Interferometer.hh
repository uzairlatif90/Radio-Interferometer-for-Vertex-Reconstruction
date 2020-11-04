#ifndef INTER_HEAD
#define INTER_HEAD

#include "IceRayTracing.cc"
#include "gsl/gsl_multimin.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "AntennaConfig.hh"
#include "TRandom3.h"

using namespace std;

namespace Interferometer{

  void GenerateChHitTimeAndCheckHits(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  bool CheckTrigger(int IgnoreCh[2][TotalAntennasRx]);
  
  void ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]);

  void FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],double ChDRTime[TotalAntennasRx]);
  
  void AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]);

  double Minimizer_f(const gsl_vector *v, void *params);

  void Minimizer(double InitialValues[3],double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx], double &FinalMinValue);
  
  int RootFinder_f(const gsl_vector * x, void *params, gsl_vector * f);
  
  void RootFinder(double InitialValues[3],double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],int &Status); 
  
}
#endif
