#include "Interferometer.hh"

void Interferometer::GenerateChHitTimeAndCheckHits(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]){

  double AntennaCoordTx[3]={TxCor[0],TxCor[1],-TxCor[2]};
  
  int DHits=0,RHits=0;
  // double LangD=0;double LangR=0;double LangRa=0;
  double timeD[TotalAntennasRx], timeR[TotalAntennasRx], timeRa[TotalAntennasRx];
  double RangD[TotalAntennasRx], RangR[TotalAntennasRx], RangRa[TotalAntennasRx];
  double RangRay[2][TotalAntennasRx];  
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);

    // LangD=RTresults[0];LangR=RTresults[1];LangRa=RTresults[2];
    timeD[iRx]=RTresults[3]*pow(10,9);
    timeR[iRx]=RTresults[4]*pow(10,9);
    timeRa[iRx]=RTresults[5]*pow(10,9);
    RangD[iRx]=RTresults[6];
    RangR[iRx]=RTresults[7];
    RangRa[iRx]=RTresults[8];
    
    timeRay[0][iRx]=timeD[iRx];
    timeRay[1][iRx]=timeR[iRx];

    RangRay[0][iRx]=RangD[iRx];
    RangRay[1][iRx]=RangR[iRx];
 
    if(RangD[iRx]==0 && RangRa[iRx]!=0){
      //cout<<"direct becomes refracted "<<endl;
      timeRay[0][iRx]=timeRa[iRx];
      RangRay[0][iRx]=RangRa[iRx];
    }
    if(RangR[iRx]==0 && RangRa[iRx]!=0){
      //cout<<"reflected becomes refracted "<<endl;
      timeRay[1][iRx]=timeRa[iRx];
      RangRay[1][iRx]=RangRa[iRx];
    }

    if(RangRay[0][iRx]==0){
      IgnoreCh[0][iRx]=0;
    }

    if(RangRay[1][iRx]==0){
      IgnoreCh[1][iRx]=0;
    }
    
    if(RangRay[0][iRx]!=0){
      DHits++;
    }

    if(RangRay[1][iRx]!=0){
      RHits++;
    }
    
    delete [] RTresults;  
  }  
  //cout<<"hits are "<<DHits<<" "<<RHits<<endl; 
}

bool Interferometer::CheckTrigger(int IgnoreCh[2][TotalAntennasRx]){

  int count[2][2]={{0,0},{0,0}};
  
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(IgnoreCh[iray][iRx]!=0 && iRx<4){
	count[iray][0]++;
      }
      if(IgnoreCh[iray][iRx]!=0 && iRx>3){
	count[iray][1]++;
      }
    }
  }

  bool RayTrigger[2]={false,false};

  for(int iray=0;iray<2;iray++){
    if(count[iray][0]+count[iray][1]>=3){
      if(count[iray][0]>=1 && count[iray][1]>=2){
	RayTrigger[iray]=true;
      } 
      if(count[iray][1]>=1 && count[iray][0]>=2){
	RayTrigger[iray]=true;
      }
    }
  }

  bool DidStationTrigger=false;

  if(RayTrigger[0]==true || RayTrigger[1]==true){
    DidStationTrigger=true;
  }

  if(RayTrigger[0]==false || RayTrigger[1]==false){
    DidStationTrigger=false;
  }

  return DidStationTrigger;
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

void Interferometer::FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],double ChDRTime[TotalAntennasRx]){

  double FirstHitTime[2];
  int FirstHitCh[2];
  vector <double> ChHitTimeB[2];

  for(int iray=0;iray<2;iray++){ 
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(IgnoreCh[iray][iRx]==1){
	ChHitTimeB[iray].push_back(ChHitTime[iray][iRx]);
      }
    }
  }

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(IgnoreCh[0][iRx]==1 && IgnoreCh[1][iRx]==1){
      ChDRTime[iRx]=ChHitTime[1][iRx]-ChHitTime[0][iRx];
    }else{
      ChDRTime[iRx]=0;
    }
  }

  if(ChHitTimeB[0].size()!=0 && ChHitTimeB[1].size()!=0){
    for(int iray=0;iray<2;iray++){   
      FirstHitTime[iray]=TMath::MinElement(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      FirstHitCh[iray]=TMath::LocMin(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      
      for(int iRx=0;iRx<TotalAntennasRx;iRx++){
	ChHitTime[iray][iRx]=ChHitTime[iray][iRx]-ChHitTimeB[0][FirstHitCh[0]];
      }      
    }
  }

  if(ChHitTimeB[0].size()!=0 && ChHitTimeB[1].size()==0){    
    FirstHitTime[0]=TMath::MinElement(ChHitTimeB[0].size(),ChHitTimeB[0].data());
    FirstHitCh[0]=TMath::LocMin(ChHitTimeB[0].size(),ChHitTimeB[0].data());
    
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChHitTime[0][iRx]=ChHitTime[0][iRx]-ChHitTimeB[0][FirstHitCh[0]];
    }
  }

  if(ChHitTimeB[0].size()==0 && ChHitTimeB[1].size()!=0){    
    FirstHitTime[1]=TMath::MinElement(ChHitTimeB[1].size(),ChHitTimeB[1].data());
    FirstHitCh[1]=TMath::LocMin(ChHitTimeB[1].size(),ChHitTimeB[1].data());
    
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChHitTime[1][iRx]=ChHitTime[1][iRx]-ChHitTimeB[1][FirstHitCh[1]];
    }
  }

  if(ChHitTimeB[0].size()==0 && ChHitTimeB[1].size()==0){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChHitTime[0][iRx]=0;
      ChHitTime[1][iRx]=0;
    }
  }
}

void Interferometer::AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]){

  ////Introduce jitter +-JitterNumber
  TRandom3 *RandNumIni = new TRandom3(0);
  for (int iRx=0; iRx<TotalAntennasRx; iRx++){
    for(int iray=0;iray<2;iray++){ 
      double RandNum = (RandNumIni->Rndm(iRx)*2*JitterNumber)-JitterNumber;
      ChHitTime[iray][iRx]=ChHitTime[iray][iRx]+RandNum;
    }
  }
}

double Interferometer::Minimizer_f(const gsl_vector *v, void *params){

  double x, y, z;
  double *p = (double *)params; 
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
  z = gsl_vector_get(v, 2);
  double AntennaCoordTx[3]={x,y,z};

  double timeRay[2][TotalAntennasRx];
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
    ChDRTime[iRx]=0;
  }
  
  Interferometer::GenerateChHitTimeAndCheckHits(AntennaCoordTx,timeRay,IgnoreCh);

  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);

  if(CheckStationTrigger==true){
    Interferometer::FindFirstHitAndNormalizeHitTime(timeRay,IgnoreCh,ChDRTime);
    
    double chi2=0;
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(p[2*TotalAntennasRx +iRx]!=0){
      	chi2+=pow(timeRay[0][iRx] - p[0+iRx],2);
      }
      // if(p[3*TotalAntennasRx +iRx]!=0){
      // 	chi2+=pow(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx],2);
      // }
      // if(p[2*TotalAntennasRx +iRx]!=0 && p[3*TotalAntennasRx +iRx]!=0){
      //   chi2+=pow(ChDRTime[iRx]- (p[1*TotalAntennasRx+iRx] - p[0+iRx]) ,2);
      // }
    }
    return chi2;
  }else{ 
    return GSL_NAN;
  }  
}

void Interferometer::Minimizer(double InitialValues[3],double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx], double &FinalMinValue){

  double ParameterArray[4*TotalAntennasRx];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
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
  gsl_vector_set (XYZVec, 0, InitialValues[0]);
  gsl_vector_set (XYZVec, 1, InitialValues[1]);
  gsl_vector_set (XYZVec, 2, InitialValues[2]);

  /* Set initial step sizes to 1 */
  XYZStepSizeVector = gsl_vector_alloc (3);
  gsl_vector_set (XYZStepSizeVector, 0, 10);
  gsl_vector_set (XYZStepSizeVector, 1, 10);
  gsl_vector_set (XYZStepSizeVector, 2, 10);
  
  /* Initialize method and iterate */
  MultiDimMinFunc.n = 3;
  MultiDimMinFunc.f = Minimizer_f;
  MultiDimMinFunc.params = ParameterArray;

  MinimizerWorkSpace = gsl_multimin_fminimizer_alloc (MinimisationType, 3);
  gsl_multimin_fminimizer_set (MinimizerWorkSpace, &MultiDimMinFunc, XYZVec, XYZStepSizeVector);
  
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
      //     printf ("converged to minimum at\n");
      //   }
      // printf ("%5d %7.3f %7.3f %7.3f f() = %7.3f size = %.3f\n",
      //         Iter,
      //         gsl_vector_get (MinimizerWorkSpace->x, 0),
      //         gsl_vector_get (MinimizerWorkSpace->x, 1),
      // 	      gsl_vector_get (MinimizerWorkSpace->x, 2),
      //         MinimizerWorkSpace->fval,Size);
     
      FinalMinValue=MinimizerWorkSpace->fval;
      FinalTxCor[0]=gsl_vector_get (MinimizerWorkSpace->x, 0);
      FinalTxCor[1]=gsl_vector_get (MinimizerWorkSpace->x, 1);
      FinalTxCor[2]=-gsl_vector_get (MinimizerWorkSpace->x, 2);      
    }
  while (Status == GSL_CONTINUE && Iter < 1000);

  gsl_vector_free(XYZVec);
  gsl_vector_free(XYZStepSizeVector);
  gsl_multimin_fminimizer_free (MinimizerWorkSpace);
}

int Interferometer::RootFinder_f(const gsl_vector * x, void *params, gsl_vector * f){
 
  const double X = gsl_vector_get (x, 0);
  const double Y = gsl_vector_get (x, 1);
  const double Z = gsl_vector_get (x, 2);
  double AntennaCoordTx[3]={X,Y,Z};
  
  double *p = (double *)params;
  
  double timeRay[2][TotalAntennasRx];
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
    ChDRTime[iRx]=0;
  }
  
  Interferometer::GenerateChHitTimeAndCheckHits(AntennaCoordTx,timeRay,IgnoreCh);
  
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);

  if(CheckStationTrigger==true){
    Interferometer::FindFirstHitAndNormalizeHitTime(timeRay,IgnoreCh,ChDRTime);
    
    double chi2a=0, chi2b=0, chi2c=0;

    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(p[2*TotalAntennasRx +iRx]!=0){
	chi2a+=pow(timeRay[0][iRx] - p[0+iRx],2);
      }
      if(p[3*TotalAntennasRx +iRx]!=0){
      	chi2b+=pow(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx],2);
      }
      // if(p[2*TotalAntennasRx +iRx]!=0){
      // 	chi2c+=chi2a;
      // }
      // if(p[3*TotalAntennasRx +iRx]!=0){
      // 	chi2c+=chi2b;
      // }
      if(p[2*TotalAntennasRx +iRx]!=0 && p[3*TotalAntennasRx +iRx]!=0){
        chi2c+=pow(ChDRTime[iRx]- (p[1*TotalAntennasRx+iRx] - p[0+iRx]) ,2);
      }
    }

    const double chi2A = chi2a;
    const double chi2B = chi2b;
    const double chi2C = chi2c;
    
    gsl_vector_set (f, 0, chi2A);
    gsl_vector_set (f, 1, chi2B);
    gsl_vector_set (f, 2, chi2C);
    
    return GSL_SUCCESS;
  }else{
    gsl_vector_set (f, 0, GSL_NAN);
    gsl_vector_set (f, 1, GSL_NAN);
    gsl_vector_set (f, 2, GSL_NAN);
    
    return GSL_ERANGE; 
  }
}

void Interferometer::RootFinder(double InitialValues[3],double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],int &Status){
  
  double ParameterArray[4*TotalAntennasRx+3];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;
  const size_t n = 3;
  gsl_multiroot_function f = {&RootFinder_f, n, &ParameterArray};

  double x_init[n];
  x_init[0] =InitialValues[0];
  x_init[1] =InitialValues[1];
  x_init[2] =InitialValues[2];
 
  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  gsl_vector_set (x, 2, x_init[2]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);
  
  // printf ("iter = %3u x = % .3f % .3f % .3f "
  //         "f(x) = % .3e % .3e % .3e\n",
  //         iter,
  //         gsl_vector_get (s->x, 0),
  //         gsl_vector_get (s->x, 1),
  // 	  gsl_vector_get (s->x, 2),
  //         gsl_vector_get (s->f, 0),
  //         gsl_vector_get (s->f, 1),
  // 	  gsl_vector_get (s->f, 2));

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      
      // printf ("iter = %3u x = % .3f % .3f % .3f "
      //     "f(x) = % .3e % .3e % .3e\n",
      //     iter,
      //     gsl_vector_get (s->x, 0),
      //     gsl_vector_get (s->x, 1),
      // 	  gsl_vector_get (s->x, 2),
      //     gsl_vector_get (s->f, 0),
      //     gsl_vector_get (s->f, 1),
      // 	  gsl_vector_get (s->f, 2));
      
      if (status)   /* check if solver is stuck */
        break;

      status =gsl_multiroot_test_residual (s->f, 1e-3);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  FinalTxCor[0]=gsl_vector_get (s->x, 0);
  FinalTxCor[1]=gsl_vector_get (s->x, 1);
  FinalTxCor[2]=-gsl_vector_get (s->x, 2);
  Status=status;
  
  //printf ("status = %s %3u \n", gsl_strerror (status), status);

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
}
