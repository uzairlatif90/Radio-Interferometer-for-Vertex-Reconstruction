#include "Interferometer.hh"

void Interferometer::XYZtoThPhR(Double_t XYZ[3],Double_t ThPhR[3]){
  ThPhR[0]=atan2(sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]),XYZ[2]);
  ThPhR[1]=atan2(XYZ[1],XYZ[0]);
  ThPhR[2]=sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2]);
}

void Interferometer::ThPhRtoXYZ(Double_t ThPhR[3],Double_t XYZ[3]){
  XYZ[0]=ThPhR[2]*sin(ThPhR[0])*cos(ThPhR[1]);
  XYZ[1]=ThPhR[2]*sin(ThPhR[0])*sin(ThPhR[1]);
  XYZ[2]=ThPhR[2]*cos(ThPhR[0]);     
}

void Interferometer::GenerateChHitTimeAndCheckHits(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]){

  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]};
  
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
      if(IgnoreCh[iray][iRx]!=0 && iRx<8){
	count[iray][0]++;
      }
      if(IgnoreCh[iray][iRx]!=0 && iRx>7){
	count[iray][1]++;
      }
    }
  }

  bool RayTrigger[2]={false,false};

  for(int iray=0;iray<2;iray++){
    if(count[iray][0]+count[iray][1]>=3){
      if(count[iray][0]>=3){
	RayTrigger[iray]=true;
      } 
      if(count[iray][1]>=3){
	RayTrigger[iray]=true;
      }
    }
  }

  bool DidStationTrigger=false;

  if(RayTrigger[0]==true || RayTrigger[1]==true){
    DidStationTrigger=true;
  }

  if(RayTrigger[0]==false && RayTrigger[1]==false){
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

  double FirstHitTime[2]={0,0};
  int FirstHitCh[2]={0,0};

  double LastHitTime[2]={0,0};
  int LastHitCh[2]={0,0};
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

  for(int iray=0;iray<2;iray++){
    if(ChHitTimeB[iray].size()!=0){
      FirstHitTime[iray]=TMath::MinElement(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      FirstHitCh[iray]=TMath::LocMin(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());

      LastHitTime[iray]=TMath::MaxElement(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      LastHitCh[iray]=TMath::LocMax(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
     }
   }
  
  for(int iray=0;iray<2;iray++){
    if(ChHitTimeB[iray].size()!=0){
      for(int iRx=0;iRx<TotalAntennasRx;iRx++){
	if(IgnoreCh[iray][iRx]==1){
	  ChHitTime[iray][iRx]=ChHitTime[iray][iRx]-FirstHitTime[0];
	}else{
	  ChHitTime[iray][iRx]=0;
	}
      }
    }else{
      for(int iRx=0;iRx<TotalAntennasRx;iRx++){
	ChHitTime[iray][iRx]=0;
      }
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

  double *p = (double *)params; 

  double ExpectedX=p[6*TotalAntennasRx];
  double ExpectedY=p[6*TotalAntennasRx+1];
  double ExpectedZ=p[6*TotalAntennasRx+2];
  double distance=sqrt( pow(ExpectedX,2)+ pow(ExpectedY,2)+ pow(ExpectedZ,2));

  double theta, phi, r;
  theta = gsl_vector_get(v, 0)*(Interferometer::pi/180.0);
  phi = gsl_vector_get(v, 1)*(Interferometer::pi/180.0);
  r = gsl_vector_get(v, 2);

  if(r>distance-500 && r<distance+500){

    Double_t ThPhR[3]={theta,phi,r};
    Double_t XYZ[3]={0,0,0};
    Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
    double AntennaCoordTx[3]={0,0,0};
    for(int ixyz=0;ixyz<3;ixyz++){
      AntennaCoordTx[ixyz]=XYZ[ixyz];
    }
    
    double timeRay[2][TotalAntennasRx];
    int IgnoreCh[2][TotalAntennasRx];
    double ChDRTime[TotalAntennasRx];

    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      for(int iray=0;iray<2;iray++){
	IgnoreCh[iray][iRx]=1;
      }
      ChDRTime[iRx]=0;
    }
  
    Double_t SumSNRD=0, SumSNRR=0,SumSNR=0;
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(p[2*TotalAntennasRx +iRx]!=0){
	SumSNRD+=p[4*TotalAntennasRx +iRx];
      }
      if(p[3*TotalAntennasRx +iRx]!=0){
	SumSNRR+=p[5*TotalAntennasRx +iRx];
      }
    }
    SumSNR=SumSNRD+SumSNRR;
  
    Interferometer::GenerateChHitTimeAndCheckHits(AntennaCoordTx,timeRay,IgnoreCh);
    bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
 
    //if(CheckStationTrigger==true){
    Interferometer::FindFirstHitAndNormalizeHitTime(timeRay,IgnoreCh,ChDRTime);
    
    double chi2=0;
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
 
      if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
      	chi2+=pow((p[4*TotalAntennasRx +iRx]*(timeRay[0][iRx] - p[0+iRx]))/SumSNR,2);
      }
      if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      	chi2+=pow((p[5*TotalAntennasRx +iRx]*(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx]))/SumSNR,2);
      }
      // if(p[2*TotalAntennasRx +iRx]!=0 && p[3*TotalAntennasRx +iRx]!=0){
      //   chi2+=pow((p[4*TotalAntennasRx +iRx]*p[5*TotalAntennasRx +iRx]*(ChDRTime[iRx]- (p[1*TotalAntennasRx+iRx] - p[0+iRx])))/(SumSNR*SumSNR) ,2);
      // }
    }
   
    return chi2;
    // }else{ 
    //   return GSL_NAN;
    // }
  }else{
    return GSL_NAN;
  }
  
}

void Interferometer::Minimizer(double InitialTxCor[3], double FinalTxCor[3], double ExpectedTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], double &FinalMinValue, int &Iterations){

  double ParameterArray[6*TotalAntennasRx+3];
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

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }
  ParameterArray[iEnt]=ExpectedTxCor[0];
  ParameterArray[iEnt+1]=ExpectedTxCor[1];
  ParameterArray[iEnt+2]=ExpectedTxCor[2];
  
  const gsl_multimin_fminimizer_type *MinimisationType = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *MinimizerWorkSpace = NULL;
  gsl_vector *XYZStepSizeVector, *XYZVec;
  gsl_multimin_function MultiDimMinFunc;

  size_t Iter = 0;
  int Status;
  double Size;
    
  /* Starting point */
  XYZVec = gsl_vector_alloc (3);
  gsl_vector_set (XYZVec, 0, InitialTxCor[0]);
  gsl_vector_set (XYZVec, 1, InitialTxCor[1]);
  gsl_vector_set (XYZVec, 2, InitialTxCor[2]);

  double ExpectedDisplacement=sqrt(pow(ExpectedTxCor[0],2) + pow(ExpectedTxCor[1],2)+ pow(ExpectedTxCor[2],2) );
  double ThetaPhiStepSize=atan(ExpectedUncertainty/ExpectedDisplacement)*(180./Interferometer::pi);
  
  /* Set initial step sizes to 1 */
  XYZStepSizeVector = gsl_vector_alloc (3);
  gsl_vector_set (XYZStepSizeVector, 0, ThetaPhiStepSize*4);
  gsl_vector_set (XYZStepSizeVector, 1, ThetaPhiStepSize*4);
  gsl_vector_set (XYZStepSizeVector, 2, ExpectedUncertainty*2);
  
  /* Initialize method and iterate */
  MultiDimMinFunc.n = 3;
  MultiDimMinFunc.f = Minimizer_f;
  MultiDimMinFunc.params = ParameterArray;

  MinimizerWorkSpace = gsl_multimin_fminimizer_alloc (MinimisationType, 3);
  gsl_multimin_fminimizer_set (MinimizerWorkSpace, &MultiDimMinFunc, XYZVec, XYZStepSizeVector);

  Iterations=0;
  do
    {
      Iter++;
      Iterations++;
      Status = gsl_multimin_fminimizer_iterate(MinimizerWorkSpace);
     
      if (Status){
        break;
      }

      Size = gsl_multimin_fminimizer_size (MinimizerWorkSpace);
      Status = gsl_multimin_test_size (Size, 1e-4);

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

      // cout<<Iter<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 0)<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 1)<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 2) <<" "<<MinimizerWorkSpace->fval<<" "<<Size<<endl;
     
      
      FinalMinValue=MinimizerWorkSpace->fval;
      FinalTxCor[0]=gsl_vector_get (MinimizerWorkSpace->x, 0);
      FinalTxCor[1]=gsl_vector_get (MinimizerWorkSpace->x, 1);
      FinalTxCor[2]=gsl_vector_get (MinimizerWorkSpace->x, 2);      
    }
  while (Status == GSL_CONTINUE && Iter < 200);
  
  gsl_vector_free(XYZVec);
  gsl_vector_free(XYZStepSizeVector);
  gsl_multimin_fminimizer_free (MinimizerWorkSpace);
}

double Interferometer::Minimizer_ftest(double InitialTxCor[3], double FinalTxCor[3], double ExpectedTxCor[3], double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],double ChSNR[2][TotalAntennasRx], double &FinalMinValue){

  double ParameterArray[6*TotalAntennasRx+3];
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

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    } 
  }
  ParameterArray[iEnt]=ExpectedTxCor[0];
  ParameterArray[iEnt+1]=ExpectedTxCor[1];
  ParameterArray[iEnt+2]=ExpectedTxCor[2];
  
  gsl_vector *XYZVec;
  /* Starting point */
  XYZVec = gsl_vector_alloc (3);
  gsl_vector_set (XYZVec, 0, InitialTxCor[0]);
  gsl_vector_set (XYZVec, 1, InitialTxCor[1]);
  gsl_vector_set (XYZVec, 2, InitialTxCor[2]);

  double chi2 =Interferometer::Minimizer_f(XYZVec, ParameterArray);
  
  gsl_vector_free(XYZVec);

  return chi2;
}

void FindApproximateMin(double RoughtInitialCor[3], double GuessResultCor[3][3], double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx],double ChSNR[2][TotalAntennasRx]){
  
  double ChDRTime[TotalAntennasRx];

  Double_t NumBinsX=10,NumBinsY=10,NumBinsZ=10;
  Double_t StartX=RoughtInitialCor[0]-500,StartY=RoughtInitialCor[1]-500,StartZ=RoughtInitialCor[2]-500;
  Double_t StopX=RoughtInitialCor[0]+500,StopY=RoughtInitialCor[1]+500,StopZ=RoughtInitialCor[2]+500;
  Double_t StepSizeX=(StopX-StartX)/NumBinsX,StepSizeY=(StopY-StartY)/NumBinsY,StepSizeZ=(StopZ-StartZ)/NumBinsZ;
 
  vector <double> RecoPar[4];
  double ParameterArray[6*TotalAntennasRx];
   
  //Interferometer::AddGaussianJitterToHitTimes(10,ChHitTime);
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime);

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
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    } 
  }

  gsl_vector *XYZVec;
  /* Starting point */
  XYZVec = gsl_vector_alloc (3);
  
  for(double i=0; i<=NumBinsX;i++){
    double Xt=i*StepSizeX+ StartX;
    for(double j=0; j<=NumBinsY;j++){	
      double Yt=j*StepSizeY + StartY;
      for(double k=0; k<=NumBinsZ;k++){
	double Zt=k*StepSizeZ + StartZ;
	
	gsl_vector_set (XYZVec, 0, Xt);
	gsl_vector_set (XYZVec, 1, Yt);
	gsl_vector_set (XYZVec, 2, Zt);
	
	double min=Interferometer::Minimizer_f(XYZVec, ParameterArray);
	if(isnan(min)==false){
	  RecoPar[0].push_back(Xt);
	  RecoPar[1].push_back(Yt);
	  RecoPar[2].push_back(Zt);
	  RecoPar[3].push_back(min);
	}
		
      }
    }
  }
  gsl_vector_free(XYZVec);
  
  int Nmin=0;
  while(Nmin<3){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];
	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=-RecoPar[2][FinalMinValueBin];
    RecoPar[3][FinalMinValueBin]=10000000000;
    
    Nmin++;
  }
  
}

void DoInterferometery(double InitialTxCor[3], double FinalTxCor[3], double ExpectedTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx],double ChSNR[2][TotalAntennasRx], double &FinalMinValue, double &Duration,int &Iterations){ 
  auto t1b = std::chrono::high_resolution_clock::now();

  Double_t XYZ[3]={InitialTxCor[0],InitialTxCor[1],InitialTxCor[2]};
  Double_t ThPhR[3]={0,0,0}; 
  Interferometer::XYZtoThPhR(XYZ,ThPhR);
  ThPhR[0]=ThPhR[0]*(180./Interferometer::pi);
  ThPhR[1]=ThPhR[1]*(180./Interferometer::pi);
  
  const int InitialValueNum=1;
  double InitialTxCor_TPR[InitialValueNum][3]={{ThPhR[0],ThPhR[1],ThPhR[2]}}; 
  
  vector <double> RecoXYZValues[3];
  vector <double> RecofMinValues;
  
  int FinalMinValueBin=0;
  double JitterNumber=5; 
  double ChDRTime[TotalAntennasRx];
  //Interferometer::AddGaussianJitterToHitTimes(JitterNumber,ChHitTime);
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime);
  
  int iInVal=0;
  while(iInVal<InitialValueNum){
    double DummyRecoCor[3];
    double DummyMin;
    Interferometer::Minimizer(InitialTxCor_TPR[iInVal],DummyRecoCor,ExpectedTxCor,ExpectedUncertainty,ChHitTime,IgnoreCh,ChSNR,DummyMin,Iterations);
    //cout<<"Attempt No. "<<iInVal<<" :: Reco Results are: Xreco="<<RecoXYZValues[iInVal][0]<<" ,Yreco="<<RecoXYZValues[iInVal][1]<<" ,Zreco="<<RecoXYZValues[iInVal][2]<<" "<<RecofMinValues[iInVal]<<endl;
 
    RecoXYZValues[0].push_back(DummyRecoCor[0]);
    RecoXYZValues[1].push_back(DummyRecoCor[1]);
    RecoXYZValues[2].push_back(DummyRecoCor[2]);
    RecofMinValues.push_back(DummyMin);
    iInVal++;
  }
  
  if(RecofMinValues.size()!=0){
    FinalMinValueBin=TMath::LocMin(RecofMinValues.size(),RecofMinValues.data());
    FinalMinValue=RecofMinValues[FinalMinValueBin];
    FinalTxCor[0]=RecoXYZValues[0][FinalMinValueBin];
    FinalTxCor[1]=RecoXYZValues[1][FinalMinValueBin];
    FinalTxCor[2]=RecoXYZValues[2][FinalMinValueBin];
  }else{
    FinalTxCor[0]=0;
    FinalTxCor[1]=0;
    FinalTxCor[2]=0;
  }

  ThPhR[0]=FinalTxCor[0];
  ThPhR[1]=FinalTxCor[1];
  ThPhR[2]=FinalTxCor[2];
  
  ThPhR[0]=ThPhR[0]*(Interferometer::pi/180.0);
  ThPhR[1]=ThPhR[1]*(Interferometer::pi/180.0);

  
  Interferometer::ThPhRtoXYZ(ThPhR,FinalTxCor);   
  
  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();

  Duration=durationb/1000;
}
