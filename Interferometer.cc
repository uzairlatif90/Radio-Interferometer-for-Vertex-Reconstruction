#include "Interferometer.hh"

void Interferometer::GenerateChHitTime(double DummyTx[3],double ChHitTime[2][TotalAntennasRx]){
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){

    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-DummyTx[0],2) + pow(AntennaCoordRx[iRx][1]-DummyTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=DummyTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);

    for(int iray=0;iray<3;iray++){ 
      if(RTresults[6+iray]==0){
	RTresults[0+iray]=0;
	RTresults[3+iray]=0;
	RTresults[6+iray]=0;
      }
    }
    
    double timeD=RTresults[3]*pow(10,9);
    double timeR=RTresults[4]*pow(10,9);
    double timeRa=RTresults[5]*pow(10,9);
    double RangD=RTresults[6];
    double RangR=RTresults[7];
    double RangRa=RTresults[8];
    
    ChHitTime[0][iRx]=timeD;
    ChHitTime[1][iRx]=timeR;
    
    if(RangD==0 && RangRa!=0){
      ChHitTime[0][iRx]=timeRa;
    }
    if(RangR==0 && RangRa!=0){
      ChHitTime[1][iRx]=timeRa;
    }
    
    delete [] RTresults;  
  }
  
}

void Interferometer::ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]){

  ////Open the file
  std::ifstream ain(filename);
  int n1=0;////variable for counting total number of data points
  std::string line;
  double dummy[5]={0,0,0};////temporary variable for storing data values from the file  

  //Check if file is open and store data
  if(ain.is_open()){
    while (true){
      ain>>dummy[0]>>dummy[1]>>dummy[2];
      if( ain.eof() ) break;     
      cout<<n1<<" "<<dummy[0]<<" "<<dummy[1]<<" "<<dummy[2]<<endl;
      ChHitTime[0][n1]= dummy[1];
      ChHitTime[1][n1]= dummy[2];
      n1++;
    }
  }
  ain.close();
}

void Interferometer::FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx]){

  double FirstHitTime[2];
  int FirstHitCh[2]; 

  for(int iray=0;iray<2;iray++){ 
    FirstHitTime[iray]=TMath::MinElement(TotalAntennasRx,ChHitTime[iray]);
    FirstHitCh[iray]=TMath::LocMin(TotalAntennasRx,ChHitTime[iray]);
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChHitTime[iray][iRx]=ChHitTime[iray][iRx]-ChHitTime[iray][FirstHitCh[iray]];
    }
  }
}

void Interferometer::AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]){

  ////Introduce jitter +-JitterNumber
  TRandom3 *RandNumIni = new TRandom3(0);
  for (int iRx=0; iRx<TotalAntennasRx; iRx++){
    for(int iray=0;iray<2;iray++){ 
      double RandNum = (RandNumIni->Rndm(iRx)*2*JitterNumber)-JitterNumber;
      ChHitTime[iRx][iray]=ChHitTime[iRx][iray]+RandNum;
    }
  }
}

double Interferometer::MinimizeXYZ(const gsl_vector *v, void *params){

  double x, y, z;
  double *p = (double *)params;  
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
  z = gsl_vector_get(v, 2);
  
  double AntennaCoordTx[3]={x,y,z};

  // double LangD=0;double LangR=0;double LangRa=0;
  double timeD[TotalAntennasRx];
  double timeR[TotalAntennasRx];
  double timeRa[TotalAntennasRx];
  double RangD=0;
  double RangR=0;
  double RangRa=0;
  double timeRay[2][TotalAntennasRx];  
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);

    for(int iray=0;iray<3;iray++){ 
      if(RTresults[6+iray]==0){
	RTresults[0+iray]=0;
	RTresults[3+iray]=0;
	RTresults[6+iray]=0;
      }
    }   

    // LangD=RTresults[0];LangR=RTresults[1];LangRa=RTresults[2];
    timeD[iRx]=RTresults[3]*pow(10,9);
    timeR[iRx]=RTresults[4]*pow(10,9);
    timeRa[iRx]=RTresults[5]*pow(10,9);
    RangD=RTresults[6];
    RangR=RTresults[7];
    RangRa=RTresults[8];
    
    timeRay[0][iRx]=timeD[iRx];
    timeRay[1][iRx]=timeR[iRx];
    
    if(RangD==0 && RangRa!=0){
      timeRay[0][iRx]=timeRa[iRx];
    }
    if(RangR==0 && RangRa!=0){
      timeRay[1][iRx]=timeRa[iRx];
    }
    delete [] RTresults;  
  }

  Interferometer::FindFirstHitAndNormalizeHitTime(timeRay);
  
  double chi2=0;
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    chi2+=pow(timeRay[0][iRx] - p[0+iRx],2) + pow(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx],2);
  }  
  
  return chi2;

}

void Interferometer::DoTheMinimization(double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx])
{

  double OneDArrayHitTime[2*TotalAntennasRx];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      OneDArrayHitTime[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  
  
  const gsl_multimin_fminimizer_type *MinimisationType = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *MinimizerWorkSpace = NULL;
  gsl_vector *XYZStepSizeVector, *XYZVec;
  gsl_multimin_function MultiDimMinFunc;

  size_t Iter = 0;
  int Status;
  double Size;

  /* Starting point */
  XYZVec = gsl_vector_alloc (3);
  gsl_vector_set (XYZVec, 0, 1);
  gsl_vector_set (XYZVec, 0, 1);
  gsl_vector_set (XYZVec, -1, 1);

  /* Set initial step sizes to 1 */
  XYZStepSizeVector = gsl_vector_alloc (3);
  gsl_vector_set_all (XYZStepSizeVector, 1);

  /* Initialize method and iterate */
  MultiDimMinFunc.n = 3;
  MultiDimMinFunc.f = MinimizeXYZ;
  MultiDimMinFunc.params = OneDArrayHitTime;

  MinimizerWorkSpace = gsl_multimin_fminimizer_alloc (MinimisationType, 3);
  gsl_multimin_fminimizer_set (MinimizerWorkSpace, &MultiDimMinFunc, XYZVec, XYZStepSizeVector);

  double FinalMinValue=0;
  do
    {
      Iter++;
      Status = gsl_multimin_fminimizer_iterate(MinimizerWorkSpace);

      if (Status)
        break;

      Size = gsl_multimin_fminimizer_size (MinimizerWorkSpace);
      Status = gsl_multimin_test_size (Size, 1e-3);

      // if (Status == GSL_SUCCESS)
      //   {
      //     //printf ("converged to minimum at\n");
      //   }

      // // printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
      // //         Iter,
      // //         gsl_vector_get (s->x, 0),
      // //         gsl_vector_get (s->x, 1),
      // //         s->fval, size);
      FinalMinValue=MinimizerWorkSpace->fval;
      FinalTxCor[0]=gsl_vector_get (MinimizerWorkSpace->x, 0);
      FinalTxCor[1]=gsl_vector_get (MinimizerWorkSpace->x, 1);
      FinalTxCor[2]=gsl_vector_get (MinimizerWorkSpace->x, 2);
    }
  while (Status == GSL_CONTINUE && Iter < 100);
  
  gsl_vector_free(XYZVec);
  gsl_vector_free(XYZStepSizeVector);
  gsl_multimin_fminimizer_free (MinimizerWorkSpace);
  
}
