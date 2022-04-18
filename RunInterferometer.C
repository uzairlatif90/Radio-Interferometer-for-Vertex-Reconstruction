#include "Interferometer.cc"

void RunInterferometer(){
  DeclareAntennaConfig();

  double DummyTxCor[3]={-180,160,-200};
  //double DummyTxCor[3]={3,4,20};
  //double DummyTxCor[3]={0.02,3,370};
  //double DummyTxCor[3]={-500,-500,-151};
  //double DummyTxCor[3]={-500,-500,-101};
  double ExpectedTimeJitter=2;// in ns
  double ExpectedPositionUncertainty=5;//in m
  
  bool RefineRecoResults=false;
  
  vector <double> ChHitTime[2]; ////Channel Hit Time
  vector <int> IgnoreCh[2];
  vector <double> ChSNR[2]; 
  ChHitTime[0].resize(TotalAntennasRx);
  ChHitTime[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  ChSNR[0].resize(TotalAntennasRx);
  ChSNR[1].resize(TotalAntennasRx);
  
  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  double FinalMinValue=0;
  int FinalMinValueBin=0;
  int Iterations;
  
  double InitialTxCor_ThPhR[3]={0,0,0};
  double InitialTxCor_XYZ[3]={0,0,0};
  double FinalTxCor_ThPhR[3]={0,0,0};
  double FinalTxCor_XYZ[3]={0,0,0};
  double TrueTxCor_ThPhR[3]={0,0,0};
  double TrueTxCor_XYZ[3]={0,0,0};
  int IsItBelowStation;
  double ArrivalDirection[3];
  
  TRandom3 *RandNumIni = new TRandom3(0);  	  
  for(int ixyz=0;ixyz<3;ixyz++){
    double RandNum = 0;// (RandNumIni->Rndm(ixyz)*2-1)*ExpectedPositionUncertainty;
    TrueTxCor_XYZ[ixyz]=DummyTxCor[ixyz]-AvgAntennaCoordRx[ixyz];
    InitialTxCor_XYZ[ixyz]=DummyTxCor[ixyz]+RandNum;
  }  
  Interferometer::XYZtoThPhR(TrueTxCor_XYZ,TrueTxCor_ThPhR);
  TrueTxCor_ThPhR[0]=TrueTxCor_ThPhR[0]*(180./Interferometer::pi);
  TrueTxCor_ThPhR[1]=TrueTxCor_ThPhR[1]*(180./Interferometer::pi); 

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
      ChSNR[iray][iRx]=10;
    }
  }
  
  if(TrueTxCor_XYZ[2]+AvgAntennaCoordRx[2]<0){
    Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
  }else{
    Interferometer::GenerateChHitTimeAndCheckHits_Air(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
  }
  Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTime);
  
  // for(int iRx=0;iRx<TotalAntennasRx;iRx++){
  //   //cout<<"ignore channels are "<<iRx<<" "<<IgnoreCh[0][iRx]<<" "<<IgnoreCh[1][iRx]<<endl;
  //   cout<<"Hit Times are "<<iRx<<" "<<ChHitTime[0][iRx]<<" "<<ChHitTime[1][iRx]<<endl;
  // }
  cout<<" X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
  cout<<" Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
  if(CheckStationTrigger==true){
    
    double A=1.78;
    double B=-0.43;
    double C=0.0132;    
    IceRayTracing::SetA(A);
    IceRayTracing::SetB(B);
    IceRayTracing::SetC(C);      
    //Check if the station was hit from above or below
    IsItBelowStation=Interferometer::IsItAboveOrBelow(ChHitTime,IgnoreCh) ;

    // double X2a=Interferometer::GetChiSquaredThPhR(TrueTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, IsItBelowStation,0,50);    
    // cout<<"chi sqaured values are "<<X2a<<endl;
      
    Interferometer::GetRecieveAngle(ChHitTime, IgnoreCh, ChSNR, 100, ArrivalDirection);

    //Set up variables for interferometry
    double GuessResultCor[3][4]; 
    
    //Start counting time it takes to guess the initial condition
    auto t1 = std::chrono::high_resolution_clock::now();

    ////Get a guess direction
    Interferometer::GetApproximateMinThPhR(GuessResultCor,ChHitTime,IgnoreCh,ChSNR,IsItBelowStation,10);
      
    //Get a guess distance
    Interferometer::GetApproximateDistance(GuessResultCor,ChHitTime,IgnoreCh,ChSNR,IsItBelowStation,10);
      
    //Calculate the time it took to guess the initial condition
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    DurationInitialCondition=duration/1000;      
       
    //Convert initial guess to radians
    InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[2]=GuessResultCor[0][2];  
      
    //Set the radial width around guess distance value in meters    
    double MinimizerRadialWidth=500;

    //do the interferometry
    Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   

    if(RefineRecoResults==true){
      cout<<"1st try Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]*(180./Interferometer::pi)<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]*(180./Interferometer::pi)<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
      
      InitialTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[2]=FinalTxCor_ThPhR[2];  
	
      Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   
    }
      
    //Calcualte the total duration
    DurationTotal=DurationInitialCondition+DurationReconstruction;
    
    //Convert Initial guess direction from theta,phi,r to x,y,z
    InitialTxCor_XYZ[0]=0;
    InitialTxCor_XYZ[1]=0;
    InitialTxCor_XYZ[2]=0;
    Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
    InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
    InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 
      
    //Convert final reconstructed direction from theta,phi,r to x,y,z
    FinalTxCor_XYZ[0]=0;
    FinalTxCor_XYZ[1]=0;
    FinalTxCor_XYZ[2]=0;
    FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
    FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180); 
    Interferometer::ThPhRtoXYZ(FinalTxCor_ThPhR, FinalTxCor_XYZ);
    FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(180./Interferometer::pi);
    FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(180./Interferometer::pi);
  
    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
    cout<<"True values are: |  Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<DurationTotal<<endl;	

  }else{

    //If no reconstruction was performed then fill the final arrays with zeroes
    ArrivalDirection[0]=0;
    ArrivalDirection[1]=0;
    ArrivalDirection[2]=0;

    FinalTxCor_XYZ[0]=0;
    FinalTxCor_XYZ[1]=0;
    FinalTxCor_XYZ[2]=0;

    FinalTxCor_ThPhR[0]=0;
    FinalTxCor_ThPhR[1]=0;
    FinalTxCor_ThPhR[2]=0;
      
    cout<<" Number of hit channels is less than 3 so no reconstruction was performed!!!"<<endl;
  }
  
}
