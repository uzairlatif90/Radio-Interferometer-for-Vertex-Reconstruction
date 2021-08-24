#include "Interferometer.cc"

void RunInterferometer(){
  DeclareAntennaConfig();

  //double DummyTxCor[3]={-100,100,-200};
  //double DummyTxCor[3]={-500,-500,-151};
  double DummyTxCor[3]={-500,-500,-101};
  double ExpectedTimeJitter=2;// in ns
  double ExpectedPositionUncertainty=sqrt(5*5+pow(ExpectedTimeJitter*0.3*1.78,2));// in m
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];
  
  double FinalMinValue=0;
  int FinalMinValueBin=0;
 
  double InitialTxCor_ThPhR[3]={0,0,0};
  double InitialTxCor_XYZ[3]={0,0,0};
  double FinalTxCor_ThPhR[3]={0,0,0};
  double FinalTxCor_XYZ[3]={0,0,0};
  double TrueTxCor_ThPhR[3]={0,0,0};
  double TrueTxCor_XYZ[3]={0,0,0};
  
  TRandom3 *RandNumIni = new TRandom3(0);  	  
  for(int ixyz=0;ixyz<3;ixyz++){
    double RandNum = (RandNumIni->Rndm(ixyz)*2-1)*ExpectedPositionUncertainty;
    TrueTxCor_XYZ[ixyz]=DummyTxCor[ixyz];
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

  //Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTime);
  Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
  
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
  if(CheckStationTrigger==true){
    double GuessResultCor[3][3];
    Interferometer::GetApproximateMinThPhR(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
    Interferometer::GetApproximateDistance(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
  
    InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[2]=GuessResultCor[0][2];
    Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
     
    double Duration=0;
    int Iterations=0;
    
    Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ExpectedPositionUncertainty, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, Duration, Iterations);
    
    double dX=TrueTxCor_XYZ[0]-FinalTxCor_XYZ[0];
    double dY=TrueTxCor_XYZ[1]-FinalTxCor_XYZ[1];
    double dZ=TrueTxCor_XYZ[2]-FinalTxCor_XYZ[2];
    
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
  
    cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<" | X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<" | Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;
    
  }////Check Station Trigger

  
}
