#include "Interferometer.cc"

void RunInterferometer(){
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  double DummyTx[3]={20,10,-10};
  double FinalTxCor[3];
  double JitterNumber=5;

  DeclareAntennaConfig();
  //Interferometer::ReadChHitTimeFromData("ChHitTimesFromData.txt",ChHitTime);
  Interferometer::GenerateChHitTime(DummyTx,ChHitTime);
  //Interferometer::AddGaussianJitterToHitTimes(JitterNumber,ChHitTime);
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime);
  Interferometer::DoTheMinimization(FinalTxCor,ChHitTime);
  
  cout<<"Final Tx Cor are: X="<<FinalTxCor[0]<<" ,Y="<<FinalTxCor[1]<<" ,Z="<<FinalTxCor[2]<<endl;
}
