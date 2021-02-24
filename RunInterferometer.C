#include "Interferometer.cc"

void RunInterferometer(){
  DeclareAntennaConfig();

  double DummyTxCor[3]={-100,100,-200};
  double ExpectedTimeJitter=5;// in ns
  double ExpectedUncertaintyPosition=sqrt(5*5+pow(ExpectedTimeJitter*0.3*1.78,2));// in m
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];  
  
  double FinalMinValue=0;
  int FinalMinValueBin=0;
  
  double InitialTxCor[3];
  double ExpectedTxCor[3];
  double FinalTxCor[3];

  TRandom3 *RandNumIni = new TRandom3(0);
	  
  if(DummyTxCor[2]>0){
    DummyTxCor[2]=-1;
  }
	  
  for(int ixyz=0;ixyz<3;ixyz++){
    double RandNum = (RandNumIni->Rndm(ixyz)*2-1)*ExpectedUncertaintyPosition;
    ExpectedTxCor[ixyz]=DummyTxCor[ixyz];
    InitialTxCor[ixyz]=DummyTxCor[ixyz]+RandNum;
  }
	  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }
	
  Interferometer::GenerateChHitTimeAndCheckHits(ExpectedTxCor,ChHitTime,IgnoreCh);
  Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTime);

  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
  if(CheckStationTrigger==true){
    double Duration=0;
    int Iterations=0;
    DoInterferometery(InitialTxCor, FinalTxCor, ExpectedTxCor, ExpectedUncertaintyPosition, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, Duration, Iterations);

    double dX=DummyTxCor[0]-FinalTxCor[0];
    double dY=DummyTxCor[1]-FinalTxCor[1];
    double dZ=DummyTxCor[2]-FinalTxCor[2];
    
    double CheckReco=0;
    CheckReco=sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2))/3;
    cout<<"Final Reco Results are: dX="<<dX<<" ,dY="<<dY<<" ,dZ="<<dZ<<" |  Xtrue="<<DummyTxCor[0]<<" ,Ytrue="<<DummyTxCor[1]<<" ,Ztrue="<<DummyTxCor[2]<<" | Xreco="<<FinalTxCor[0]<<" ,Yreco="<<FinalTxCor[1]<<" ,Zreco="<<FinalTxCor[2]<<" | Fn Min Value: "<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;
    if(CheckReco>ExpectedUncertaintyPosition){
      cout<<"Minimization Failed!"<<endl;
    }	 

  }////Check Station Trigger

  
}
