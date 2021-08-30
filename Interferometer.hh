#ifndef INTER_HEAD
#define INTER_HEAD

#include "IceRayTracing.cc"
#include "gsl/gsl_multimin.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "AntennaConfig.hh"
#include <iostream>     // std::cout
#include <algorithm>    // std::max

using namespace std;

namespace Interferometer{

  double pi=4*atan(1);

  void XYZtoThPhR(Double_t XYZ[3],Double_t ThPhR[3]);
  
  void ThPhRtoXYZ(Double_t ThPhR[3],Double_t XYZ[3]);
  
  void GenerateChHitTimeAndCheckHits(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  void GenerateChHitTimeAndCheckHits_Cnz(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  void GenerateChHitTimeAndCheckHits_Air(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);
  
  bool CheckTrigger(int IgnoreCh[2][TotalAntennasRx]);
  
  void ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]);

  void FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],double ChDRTime[TotalAntennasRx], int ChHitOrder[TotalAntennasRx]);
  
  void AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]);

  double Minimizer_f(const gsl_vector *v, void *params);

  double Minimizer_fCnz(const gsl_vector *v, void *params);

  double MinimizerThPh(double x, void * params);
  
  void MinimizerThPhR(double InitialTxCor_XYZ[3], double InitialTxCor_ThPhR[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int ChHitOrder[TotalAntennasRx], double &FinalMinValue, int &Iterations);

  void RootThPhR(double InitialTxCor_XYZ[3], double InitialTxCor_ThPhR[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int ChHitOrder[TotalAntennasRx], double &FinalMinValue, int &Iterations);

  void SearchApproxiMin(int C_nz, double StartCor[3],double GuessResultCor[3][3],double ParameterArray[6*TotalAntennasRx+12],int &iEnt);

  void GetApproximateMinUserCor(vector <double> UserCor[3] ,double GuessResultCor[3][3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx]);
  
  void GetApproximateMinThPhR(double GuessResultCor[3][3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx]);

  void GetApproximateDistance(double GuessResultCor[3][3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx]);
  
  void DoInterferometery(double InitialTxCor[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx],double ChSNR[2][TotalAntennasRx], double &FinalMinValue, double &Duration,int &Iterations);

}
#endif
