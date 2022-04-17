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

  //expected uncertainty in a verttex position which we are trying to reconstruct. Leave this unchanged...
  double ExpectedUncertainty=5;//in m

  void XYZtoThPhR(Double_t XYZ[3],Double_t ThPhR[3]);
  
  void ThPhRtoXYZ(Double_t ThPhR[3],Double_t XYZ[3]);
  
  void GenerateChHitTimeAndCheckHits(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  void GenerateChHitTimeAndCheckHits_Cnz(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  void GenerateChHitTimeAndCheckHits_Air(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);
  
  bool CheckTrigger(int IgnoreCh[2][TotalAntennasRx]);
  
  void ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]);

  void FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],int NormalizeRay);
  
  void AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]);

  int IsItAboveOrBelow(double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]);

  double Minimizer_f1D(double m, void *params);

  double Minimizer_f(const gsl_vector *v, void *params);

  double Minimizer_fCnz(const gsl_vector *v, void *params);

  double MinimizerThPh(double x, void * params);

  void MinimizerThPhR1D(double m, double a, double b, double &FinalMinValue, double &FnMinValue, int &Iterations, void *parameters);
  
  void MinimizerThPhR(double InitialTxCor_XYZ[3], double InitialTxCor_ThPhR[3], double FinalTxCor[3], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], double &FinalMinValue, int &Iterations, double MinimizerRadialWidth,int IsItBelowStation, int max_iter);

  double GetChiSquaredThPhR(double UserCor[3], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int IsItBelowStation, int CnzOrEnz, int max_iter);

  void GetRecieveAngle( double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int max_iter, double StartCor[3]);

  void SearchApproxiMin(double GuessResultCor[3][4],double ParameterArray[6*TotalAntennasRx+14], bool CheckAboveSurface);

  void GetApproximateMinUserCor(vector <double> UserCor[3] ,double GuessResultCor[3][4], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int IsItBelowStation, int max_iter);
  
  void GetApproximateMinThPhR(double GuessResultCor[3][4], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int IsItBelowStation, int max_iter);

  void GetApproximateDistance(double GuessResultCor[3][4], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int IsItBelowStation, int max_iter);

  void GetRecoFixedR(double InitialTxCor[3], double FinalTxCor[3], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], double &FixedR, int IsItBelowStation, int max_iter);
  
  void DoInterferometery(double InitialTxCor[3], double FinalTxCor[3], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx],double ChSNR[2][TotalAntennasRx], double &FinalMinValue, double &Duration,int &Iterations, double MinimizerRadialWidth, int IsItBelowStation, int max_iter);

}
#endif
