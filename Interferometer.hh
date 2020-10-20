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

  bool GenerateChHitTimeAndCheckHits(double DummyTx[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  bool CheckTrigger(int IgnoreCh[2][TotalAntennasRx]);
  
  void ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]);

  void FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx]);

  void AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]);

  double MinimizeXYZ(const gsl_vector *v, void *params);

  void DoTheMinimization(double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx], double &FinalMinValue);
  
  int TxXYZRoots(const gsl_vector * x, void *params, gsl_vector * f);
  
  void print_state (size_t iter, gsl_multiroot_fsolver * s);
  
  void FindRoots3D(double TxCor[3],double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx], int AdjustIniCond,int &Status);

  
}
#endif
