#include "../Interferometer/Interferometer.cc"

void ReconstructSIM(int eventID, double Theta, double Phi, double R){
  DeclareAntennaConfigARA(2);

  if(Phi>=180){
    Phi=Phi-360;
  }
  double ExpectedTimeJitter=5;// in ns
  double ExpectedPositionUncertainty=5;// in m
  
  double ChDRTime[TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];  
  
  int FinalMinValueBin=0; 
   
  double eventNum;
  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  double FinalMinValue;
  int Iterations;
  double FinalTxCor_XYZ[3];
  double FinalTxCor_ThPhR[3];
  double FinalTxCor_XYZ_fR[3];
  double FinalTxCor_ThPhR_fR[3];
  double InitialTxCor_XYZ[3];
  double InitialTxCor_ThPhR[3];
  double TrueTxCor_XYZ[3];
  double TrueTxCor_ThPhR[3];
  
  double ChHitTime[2][16];
  int IgnoreCh[2][16];
  double dXYZ[3];
  double dThPhR[3];

  TString OutputFileName="./output/";
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
  RecoTree->Branch("FinalTxCor_XYZ_fR",FinalTxCor_XYZ_fR,"FinalTxCor_XYZ_fR[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR_fR",FinalTxCor_ThPhR_fR,"FinalTxCor_ThPhR_fR[3]/D");
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
    double RandNum = (RandNumIni->Rndm(ixyz)*2-1)*ExpectedPositionUncertainty;
    InitialTxCor_XYZ[ixyz]=TrueTxCor_XYZ[ixyz]+RandNum;
  }
 	  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
      ChSNR[iray][iRx]=10;
    }
  }
  
  Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
  Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTime);
  
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
  if(CheckStationTrigger==true){
    cout<<"XYZ Tx coordinates are  "<<TrueTxCor_XYZ[0]<<" "<<TrueTxCor_XYZ[1]<<" "<<TrueTxCor_XYZ[2]<<endl;
    double GuessResultCor[3][3];
    double StartDistance=0;
    auto t1 = std::chrono::high_resolution_clock::now();

    Interferometer::GetApproximateMinThPhR(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR,StartDistance);
    Interferometer::GetApproximateDistance(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);

    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    DurationInitialCondition=duration/1000;
  
    InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[2]=GuessResultCor[0][2];
    Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
	    
    double Duration=0;
    int Iterations=0;  
    double MinimizerRadialWidth=100;
    Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ExpectedPositionUncertainty, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth);
    
    DurationTotal=DurationInitialCondition+DurationReconstruction;

    Interferometer::XYZtoThPhR(InitialTxCor_XYZ,InitialTxCor_ThPhR);
    InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
    InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 

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
    
    cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<" | X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<" | Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<DurationTotal<<endl;
  }else{

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
