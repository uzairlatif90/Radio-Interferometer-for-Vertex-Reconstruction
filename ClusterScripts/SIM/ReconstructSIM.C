#include "/data/user/ulatif/Interferometer/Interferometer.cc"

void ReconstructSIM(int eventID, double Theta, double Phi, double R){
  DeclareAntennaConfigARA(2);

  TString OutputFileName="/data/user/ulatif/SIM_Inter/output/";
  OutputFileName+="RunSim";
  OutputFileName+="Event";
  OutputFileName+=eventID;
  OutputFileName+="Theta";
  OutputFileName+=(int)Theta;
  OutputFileName+="Phi";
  OutputFileName+=(int)Phi;
  OutputFileName+="R";
  OutputFileName+=(int)R;
  OutputFileName+=".root";

  if(Phi>180){
    Phi=Phi-360;
  }
  
  double ExpectedTimeJitter=5;// in ns
  double ExpectedPositionUncertainty=5;// in m
  bool RefineRecoResults=false;

  double ChSNR[2][TotalAntennasRx];  
   
  double eventNum;
  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  double FinalMinValue;
  int FinalMinValueBin=0; 
  int Iterations;
  double FinalTxCor_XYZ[3];
  double FinalTxCor_ThPhR[3];
  double InitialTxCor_XYZ[3];
  double InitialTxCor_ThPhR[3];
  double TrueTxCor_XYZ[3];
  double TrueTxCor_ThPhR[3];

  int IsItBelowStation;
  double ArrivalDirection[3];
 
  
  double ChHitTime[2][16];
  int IgnoreCh[2][16];
  double dXYZ[3];
  double dThPhR[3];
 
  ///Create the output root file
  TFile *OutputFile=new TFile(OutputFileName,"RECREATE"); 
  
  TTree *RecoTree = new TTree("RecoTree","Reco info about Event");
  RecoTree->Branch("eventNum",&eventNum,"eventNum/D");
  RecoTree->Branch("DurationTotal",&DurationTotal,"DurationTotal/D");
  RecoTree->Branch("DurationReconstruction",&DurationReconstruction,"DurationReconstruction/D");
  RecoTree->Branch("DurationInitialCondition",&DurationInitialCondition,"DurationInitialCondition/D");
  RecoTree->Branch("FinalMinValue",&FinalMinValue,"FinalMinValue/D");
  RecoTree->Branch("Iterations",&Iterations,"Iterations/I");
 
  RecoTree->Branch("FinalTxCor_XYZ",FinalTxCor_XYZ,"FinalTxCor_XYZ[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR",FinalTxCor_ThPhR,"FinalTxCor_ThPhR[3]/D");
  RecoTree->Branch("InitialTxCor_XYZ",InitialTxCor_XYZ,"InitialTxCor_XYZ[3]/D");
  RecoTree->Branch("InitialTxCor_ThPhR",InitialTxCor_ThPhR,"InitialTxCor_ThPhR[3]/D");  
  RecoTree->Branch("TrueTxCor_XYZ",TrueTxCor_XYZ,"TrueTxCor_XYZ[3]/D");
  RecoTree->Branch("TrueTxCor_ThPhR",TrueTxCor_ThPhR,"TrueTxCor_ThPhR[3]/D");  
  
  RecoTree->Branch("ChHitTime",ChHitTime,"ChHitTime[2][16]/D");
  RecoTree->Branch("IgnoreCh",IgnoreCh,"IgnoreCh[2][16]/I");  
  RecoTree->Branch("dXYZ",dXYZ,"dXYZ[3]/D");
  RecoTree->Branch("dThPhR",dThPhR,"dThPhR[3]/D");

  eventNum=eventID;
  
  TrueTxCor_ThPhR[0]=Theta*(Interferometer::pi/180.0);
  TrueTxCor_ThPhR[1]=Phi*(Interferometer::pi/180.0);
  TrueTxCor_ThPhR[2]=R;
  TrueTxCor_XYZ[0]=0;
  TrueTxCor_XYZ[1]=0;
  TrueTxCor_XYZ[2]=0;
 
  Interferometer::ThPhRtoXYZ(TrueTxCor_ThPhR,TrueTxCor_XYZ);
  TrueTxCor_ThPhR[0]=Theta;
  TrueTxCor_ThPhR[1]=Phi;
  
  TRandom3 *RandNumIni = new TRandom3(0); 
  for(int ixyz=0;ixyz<3;ixyz++){
    double RandNum = 0;//(RandNumIni->Rndm(ixyz)*2-1)*ExpectedPositionUncertainty;
    InitialTxCor_XYZ[ixyz]=TrueTxCor_XYZ[ixyz]+RandNum;
  }
 	  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){

      IgnoreCh[iray][iRx]=1;
      ChSNR[iray][iRx]=10;
    }
    // if(Theta<45){
    //   IgnoreCh[1][iRx]=0;
    // }

  }

  vector <double> ChHitTimev[2]; ////Channel Hit Time
  vector <int> IgnoreChv[2];
  vector <double> ChSNRv[2]; 
  ChHitTimev[0].resize(TotalAntennasRx);
  ChHitTimev[1].resize(TotalAntennasRx);
  IgnoreChv[0].resize(TotalAntennasRx);
  IgnoreChv[1].resize(TotalAntennasRx);
  ChSNRv[0].resize(TotalAntennasRx);
  ChSNRv[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      ChHitTimev[iray][iRx]=ChHitTime[iray][iRx];
      IgnoreChv[iray][iRx]=IgnoreCh[iray][iRx];
      ChSNRv[iray][iRx]=ChSNR[iray][iRx];
    }
  }

  if(TrueTxCor_XYZ[2]+AvgAntennaCoordRx[2]<0){
    Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTimev,IgnoreChv);  
  }else{
    Interferometer::GenerateChHitTimeAndCheckHits_Air(TrueTxCor_XYZ,ChHitTimev,IgnoreChv);  
  }  
  Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTimev);  
  
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreChv);
  if(CheckStationTrigger==true){
    cout<<"XYZ Tx coordinates are  "<<TrueTxCor_XYZ[0]<<" "<<TrueTxCor_XYZ[1]<<" "<<TrueTxCor_XYZ[2]<<endl;
        
    double A=1.78;
    double B=-0.43;
    double C=0.0132;    
    IceRayTracing::SetA(A);
    IceRayTracing::SetB(B);
    IceRayTracing::SetC(C);      
    //Check if the station was hit from above or below
    IsItBelowStation=Interferometer::IsItAboveOrBelow(ChHitTimev,IgnoreChv) ;

    // double X2a=Interferometer::GetChiSquaredThPhR(TrueTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, IsItBelowStation,0,50);    
    // cout<<"chi sqaured values are "<<X2a<<endl;
      
    Interferometer::GetRecieveAngle(ChHitTimev, IgnoreChv, ChSNRv, 100, ArrivalDirection);

    //Set up variables for interferometry
    double GuessResultCor[3][4]; 
    
    //Start counting time it takes to guess the initial condition
    auto t1 = std::chrono::high_resolution_clock::now();

    ////Get a guess direction
    Interferometer::GetApproximateMinThPhR(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,10);
      
    //Get a guess distance
    Interferometer::GetApproximateDistance(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,10);
      
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
    Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTimev, IgnoreChv, ChSNRv, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   

    if(RefineRecoResults==true){
      cout<<"1st try Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]*(180./Interferometer::pi)<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]*(180./Interferometer::pi)<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
      
      InitialTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[2]=FinalTxCor_ThPhR[2];  
	
      Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTimev, IgnoreChv, ChSNRv, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   
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

    for(int i=0;i<3;i++){
      dXYZ[i]=InitialTxCor_XYZ[i]-FinalTxCor_XYZ[i]; 
      dThPhR[i]=InitialTxCor_ThPhR[i]-FinalTxCor_ThPhR[i]; 
    }
    
    //cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<" | X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<" | Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
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

    Interferometer::XYZtoThPhR(InitialTxCor_XYZ,InitialTxCor_ThPhR);
    InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
    InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 
    
    cout<<" Initial condition is in shadow zone no reconstruction was performed!!!"<<endl;
    
  }

  RecoTree->Fill();
  RecoTree->Write();

  OutputFile->Write();
  OutputFile->Close();
  
}
