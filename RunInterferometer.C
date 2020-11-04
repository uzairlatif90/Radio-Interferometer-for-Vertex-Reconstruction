#include "Interferometer.cc"

void RunInterferometer(){
  DeclareAntennaConfig();
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx]; ////Channel Hit Time
  double ChDRTime[TotalAntennasRx]; ////Channel Hit Time
  
  double FinalTxCor[3];
  double JitterNumber=5;
 
  double FinalMinValue;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;

  double i=-10;////Xtrue
  double j=20;////Ytrue
  double k=-10;////Ztrue
  
  double DummyTx[3]={i,j,-k};
	  
  double InitialValues[3]={1,1,1};
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<3;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }
	
  Interferometer::GenerateChHitTimeAndCheckHits(DummyTx,ChHitTime,IgnoreCh);
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
  if(CheckStationTrigger==true){
    int status;
    //Interferometer::ReadChHitTimeFromData("ChHitTimesFromData.txt",ChHitTime);
    //Interferometer::AddGaussianJitterToHitTimes(JitterNumber,ChHitTime);
    Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime);
	    
    Interferometer::Minimizer(InitialValues,FinalTxCor,ChHitTime,IgnoreCh, FinalMinValue);
    cout<<"First Attempt Tx Cor are: dX="<<i-FinalTxCor[0]<<" ,dY="<<j-FinalTxCor[1]<<" ,dZ="<<k-FinalTxCor[2]<<" |  Xtrue="<<i<<" ,Ytrue="<<j<<" ,Ztrue="<<k<<" | Xreco="<<FinalTxCor[0]<<" ,Yreco="<<FinalTxCor[1]<<" ,Zreco="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
    double CheckReco=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
    if(CheckReco>=3){
      InitialValues[0]=0;
      InitialValues[1]=0;
      Interferometer::Minimizer(InitialValues,FinalTxCor,ChHitTime,IgnoreCh,FinalMinValue);
      CheckReco=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
      cout<<"Second Attempt Tx Cor are: dX="<<i-FinalTxCor[0]<<" ,dY="<<j-FinalTxCor[1]<<" ,dZ="<<k-FinalTxCor[2]<<" |  Xtrue="<<i<<" ,Ytrue="<<j<<" ,Ztrue="<<k<<" | Xreco="<<FinalTxCor[0]<<" ,Yreco="<<FinalTxCor[1]<<" ,Zreco="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
      if(CheckReco>=3){
	cout<<"Minimization Failed!"<<endl;
	countRecoFailed++;
      }
    }
	    
    CheckReco=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
    if(CheckReco<3){
      countRecoSuccess++;
    }
	 
    countTotal++;
  }////Check Station Trigger

  cout<<"SuccessReco: "<<countRecoSuccess<<", FailedReco: "<<countRecoFailed<<", TotalReco: "<<countTotal<<endl;
  
}
