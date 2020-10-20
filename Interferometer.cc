#include "Interferometer.hh"

bool Interferometer::GenerateChHitTimeAndCheckHits(double DummyTx[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]){

  double RangD[TotalAntennasRx];
  double RangR[TotalAntennasRx];
  double RangRa[TotalAntennasRx];
  double RangRay[2][TotalAntennasRx];
  bool NoChHit=true;
  int DHits=0;
  int RHits=0;
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){

    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-DummyTx[0],2) + pow(AntennaCoordRx[iRx][1]-DummyTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=DummyTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);

    for(int iray=0;iray<3;iray++){ 
      // if(RTresults[6+iray]==0){
      // 	RTresults[0+iray]=0;
      // 	RTresults[3+iray]=0;
      // 	RTresults[6+iray]=0;
      // }
    }
    
    double timeD=RTresults[3]*pow(10,9);
    double timeR=RTresults[4]*pow(10,9);
    double timeRa=RTresults[5]*pow(10,9);
    RangD[iRx]=RTresults[6];
    RangR[iRx]=RTresults[7];
    RangRa[iRx]=RTresults[8];  

    ChHitTime[0][iRx]=timeD;
    ChHitTime[1][iRx]=timeR;

    RangRay[0][iRx]=RangD[iRx];
    RangRay[1][iRx]=RangR[iRx];

    //cout<<iRx<<" "<<timeD<<" "<<timeR<<" "<<timeRa<<endl;
    
    if(RangD[iRx]==0 && RangRa[iRx]!=0){
      //cout<<"direct becomes refracted "<<endl;
      ChHitTime[0][iRx]=timeRa;
      RangRay[0][iRx]=RangRa[iRx];
    }
    if(RangR[iRx]==0 && RangRa[iRx]!=0){
      //cout<<"reflected becomes refracted "<<endl;
      ChHitTime[1][iRx]=timeRa;
      RangRay[1][iRx]=RangRa[iRx];
    }

    if(RangRay[0][iRx]==0){
      //ChHitTime[0][iRx]=100000000;
      IgnoreCh[0][iRx]=0;
    }

    if(RangRay[1][iRx]==0){
      //ChHitTime[1][iRx]=100000000;
      IgnoreCh[1][iRx]=0;
    }   
    
    if(RangRay[0][iRx]!=0 || RangRay[1][iRx]!=0){
      NoChHit=false;
    }

    if(RangRay[0][iRx]!=0){
      DHits++;
    }

    if(RangRay[1][iRx]!=0){
      RHits++;
    }
    
    delete [] RTresults;  
  }
  // if(DHits<=5 && DHits>0){
  //   cout<<"hits are "<<DHits<<" "<<RHits<<endl;
  // }
  return NoChHit;
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

void Interferometer::FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx]){

  double FirstHitTime[2];
  int FirstHitCh[2]; 

  for(int iray=0;iray<2;iray++){ 
    FirstHitTime[iray]=TMath::MinElement(TotalAntennasRx,ChHitTime[iray]);
    FirstHitCh[iray]=TMath::LocMin(TotalAntennasRx,ChHitTime[iray]);
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChHitTime[iray][iRx]=ChHitTime[iray][iRx]-ChHitTime[0][FirstHitCh[0]];
      //cout<<iray<<" "<<iRx<<" "<<ChHitTime[iray][iRx]<<endl;
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

double Interferometer::MinimizeXYZ(const gsl_vector *v, void *params){

  double x, y, z;
  double *p = (double *)params;  
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
  z = gsl_vector_get(v, 2);
  
  double AntennaCoordTx[3]={x,y,-z};

  // double LangD=0;double LangR=0;double LangRa=0;
  double timeD[TotalAntennasRx];
  double timeR[TotalAntennasRx];
  double timeRa[TotalAntennasRx];
  double RangD[TotalAntennasRx];
  double RangR[TotalAntennasRx];
  double RangRa[TotalAntennasRx];

  double timeRay[2][TotalAntennasRx];
  double RangRay[2][TotalAntennasRx];  
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);

    // for(int iray=0;iray<3;iray++){ 
    //   if(RTresults[6+iray]==0){
    // 	// RTresults[0+iray]=0;
    // 	// RTresults[3+iray]=0;
    // 	// RTresults[6+iray]=0;
    //   }
    // }   

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
      //cout<<"direct becomes refracted "<<endl;
      timeRay[1][iRx]=timeRa[iRx];
      RangRay[1][iRx]=RangRa[iRx];
    }
    delete [] RTresults;  
  }

  Interferometer::FindFirstHitAndNormalizeHitTime(timeRay);
  
  double chi2=0;
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    // if(p[2*TotalAntennasRx +iRx]!=0){
    //   chi2+=pow(timeRay[0][iRx] - p[0+iRx],2) + pow(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx],2);
    // }
    if(p[2*TotalAntennasRx +iRx]!=0){
      chi2+=pow(timeRay[0][iRx] - p[0+iRx],2);
    }
    if(p[3*TotalAntennasRx +iRx]!=0){
      chi2+=pow(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx],2);
    }
  }

  // if(chi2==0){
  //   cout<<"it is zero "<<endl;
  //   chi2=10000000;
  // }
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    //cout<<p[3*TotalAntennasRx +iRx]<<" "<<timeRay[1][iRx]<<" "<<p[1*TotalAntennasRx+iRx]<<endl;;
  }
  return chi2;

}

void Interferometer::DoTheMinimization(double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx], double &FinalMinValue){

  double OneDArrayHitTime[4*TotalAntennasRx];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      OneDArrayHitTime[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      OneDArrayHitTime[iEnt]=IgnoreCh[iray][iRx];
      //cout<<iray<<" "<<iRx<<" "<<IgnoreCh[iray][iRx]<<endl;
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
  gsl_vector_set (XYZVec, 1, 1);
  gsl_vector_set (XYZVec, 2, 1);

  /* Set initial step sizes to 1 */
  XYZStepSizeVector = gsl_vector_alloc (3);
  //gsl_vector_set_all (XYZStepSizeVector, 10);
  gsl_vector_set (XYZStepSizeVector, 0, 10);
  gsl_vector_set (XYZStepSizeVector, 1, 10);
  gsl_vector_set (XYZStepSizeVector, 2, 10);
  
  /* Initialize method and iterate */
  MultiDimMinFunc.n = 3;
  MultiDimMinFunc.f = MinimizeXYZ;
  MultiDimMinFunc.params = OneDArrayHitTime;

  MinimizerWorkSpace = gsl_multimin_fminimizer_alloc (MinimisationType, 3);
  gsl_multimin_fminimizer_set (MinimizerWorkSpace, &MultiDimMinFunc, XYZVec, XYZStepSizeVector);
  
  do
    {
      Iter++;
      Status = gsl_multimin_fminimizer_iterate(MinimizerWorkSpace);

      if (Status)
        break;

      Size = gsl_multimin_fminimizer_size (MinimizerWorkSpace);
      Status = gsl_multimin_test_size (Size, 1e-2);

      if (Status == GSL_SUCCESS)
        {
          //printf ("converged to minimum at\n");
        }

      // printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
      //         Iter,
      //         gsl_vector_get (s->x, 0),
      //         gsl_vector_get (s->x, 1),
      //         s->fval, size);
      FinalMinValue=MinimizerWorkSpace->fval;
      FinalTxCor[0]=gsl_vector_get (MinimizerWorkSpace->x, 0);
      FinalTxCor[1]=gsl_vector_get (MinimizerWorkSpace->x, 1);
      FinalTxCor[2]=-gsl_vector_get (MinimizerWorkSpace->x, 2);

      //cout<<FinalMinValue<<" "<<FinalTxCor[0]<<" "<<FinalTxCor[1]<<" "<<FinalTxCor[2]<<endl;
    }
  while (Status == GSL_CONTINUE && Iter < 1000);
  
  gsl_vector_free(XYZVec);
  gsl_vector_free(XYZStepSizeVector);
  gsl_multimin_fminimizer_free (MinimizerWorkSpace);
  
}


int Interferometer::TxXYZRoots(const gsl_vector * x, void *params, gsl_vector * f){
 
  const double X = gsl_vector_get (x, 0);
  const double Y = gsl_vector_get (x, 1);
  const double Z = gsl_vector_get (x, 2);
  double *p = (double *)params;
  
  double AntennaCoordTx[3]={X,Y,Z};

  // double LangD=0;double LangR=0;double LangRa=0;
  double timeD[TotalAntennasRx];
  double timeR[TotalAntennasRx];
  double timeRa[TotalAntennasRx];
  double RangD[TotalAntennasRx];
  double RangR[TotalAntennasRx];
  double RangRa[TotalAntennasRx];

  double timeRay[2][TotalAntennasRx];
  double RangRay[2][TotalAntennasRx];  
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);

    // for(int iray=0;iray<3;iray++){ 
    //   if(RTresults[6+iray]==0){
    // 	// RTresults[0+iray]=0;
    // 	// RTresults[3+iray]=0;
    // 	// RTresults[6+iray]=0;
    //   }
    // }   

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
      //cout<<"direct becomes refracted "<<endl;
      timeRay[1][iRx]=timeRa[iRx];
      RangRay[1][iRx]=RangRa[iRx];
    }
    delete [] RTresults;  
  }

  Interferometer::FindFirstHitAndNormalizeHitTime(timeRay);
  
  double chi2a=0, chi2b=0, chi2c=0;
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(p[2*TotalAntennasRx +iRx]!=0){
      chi2a+=pow(timeRay[0][iRx] - p[0+iRx],2);
      //chi2a+=abs(timeRay[0][iRx] - p[0+iRx]);
    }
    if(p[3*TotalAntennasRx +iRx]!=0){
      chi2b+=pow(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx],2);
      //chi2b+=abs(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx]);
    }
    if(p[2*TotalAntennasRx +iRx]!=0){
      chi2c+=chi2a;
    }
    if(p[3*TotalAntennasRx +iRx]!=0){
      chi2c+=chi2b;
    }
    
  }

  // if(chi2==0){
  //   cout<<"it is zero "<<endl;
  //   chi2=10000000;
  // }
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    //cout<<p[3*TotalAntennasRx +iRx]<<" "<<timeRay[1][iRx]<<" "<<p[1*TotalAntennasRx+iRx]<<endl;;
  }

  const double chi2A = chi2a;
  const double chi2B = chi2b;
  const double chi2C = chi2c;

  gsl_vector_set (f, 0, chi2A);
  gsl_vector_set (f, 1, chi2B);
  gsl_vector_set (f, 2, chi2C);

  return GSL_SUCCESS;
}
  
void Interferometer::print_state (size_t iter, gsl_multiroot_fsolver * s){
  printf ("iter = %3u x = % .3f % .3f % .3f "
          "f(x) = % .3e % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
	  gsl_vector_get (s->f, 2));
 }

void Interferometer::FindRoots3D(double TxCor[3],double FinalTxCor[3],double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx], int AdjustIniCond, int &Status){
  
  double OneDArrayHitTime[4*TotalAntennasRx+3];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      OneDArrayHitTime[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      OneDArrayHitTime[iEnt]=IgnoreCh[iray][iRx];
      //cout<<iray<<" "<<iRx<<" "<<IgnoreCh[iray][iRx]<<endl;
      iEnt++;
    } 
  }

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;
  const size_t n = 3;
  gsl_multiroot_function f = {&TxXYZRoots, n, &OneDArrayHitTime};

  double x_init[n];

  if(AdjustIniCond==1){
    x_init[0] =TxCor[0];
    x_init[1] =TxCor[1];
    x_init[2] =TxCor[2];
  }
  if(AdjustIniCond==0){
    x_init[0] =1;
    x_init[1] =1;
    x_init[2] =-1;
  }
  //double x_init[n] = {1, 1, -1};

  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);
  gsl_vector_set (x, 2, x_init[2]);

  T = gsl_multiroot_fsolver_hybrid;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  //print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      //print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status =gsl_multiroot_test_residual (s->f, 1e-2);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  FinalTxCor[0]=gsl_vector_get (s->x, 0);
  FinalTxCor[1]=gsl_vector_get (s->x, 1);
  FinalTxCor[2]=gsl_vector_get (s->x, 2);
  
  //printf ("status = %s %3u \n", gsl_strerror (status), status);

  Status=status;
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return 0;
}
