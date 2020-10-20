#include "Interferometer.cc"

void RunInterferometer(){
  DeclareAntennaConfig();
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx]; ////Channel Hit Time

  double FinalTxCor[3];
  double JitterNumber=2;

  double FinalMinValue;
  double DummyTx[3]={10,-20,-10};
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<3;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }
  
  Interferometer::GenerateChHitTimeAndCheckHits(DummyTx,ChHitTime,IgnoreCh);
  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
  if(CheckStationTrigger==true){
    int status;
    //Interferometer::{ReadChHitTimeFromData("ChHitTimesFromData.txt",ChHitTime);
    //Interferometer::AddGaussianJitterToHitTimes(JitterNumber,ChHitTime);
    Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime);
    
    Interferometer::DoTheMinimization(FinalTxCor,ChHitTime,IgnoreCh, FinalMinValue);
    cout<<"First Tx Cor are: dX="<<DummyTx[0]-FinalTxCor[0]<<" ,dY="<<DummyTx[1]-FinalTxCor[1]<<" ,dZ="<<DummyTx[2]-FinalTxCor[2]<<"  Xi="<<DummyTx[0]<<" ,Yi="<<DummyTx[1]<<" ,Zi="<<DummyTx[2]<<"   Xi="<<FinalTxCor[0]<<" ,Yi="<<FinalTxCor[1]<<" ,Zi="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
    double checkdiff=fabs(DummyTx[0]-FinalTxCor[0])+fabs(DummyTx[1]-FinalTxCor[1])+fabs(DummyTx[2]-FinalTxCor[2]);
    if(checkdiff>3){
      Interferometer::FindRoots3D(DummyTx,FinalTxCor,ChHitTime,IgnoreCh,0,status);
     
      if(status!=0){
	Interferometer::FindRoots3D(FinalTxCor,FinalTxCor,ChHitTime,IgnoreCh,1,status);
      }
      
      if(status!=0){
	for(double itest=-500; itest<501;itest=itest+90){
	  for(double jtest=-500; jtest<501;jtest=jtest+90){
	    if(itest!=0 || jtest!=0){
	      for(double ktest=-1; ktest>-1001;ktest=ktest-40){
		
		double DummyTxtest[3]={itest,jtest,ktest};
		//cout<<"Tx Cor are: i="<<itest<<" ,j="<<jtest<<" ,k="<<ktest<<" ,status "<<status<<endl;
		Interferometer::FindRoots3D(DummyTxtest,FinalTxCor,ChHitTime,IgnoreCh,1,status);
		if(status==0){
		  itest=1000;
		  jtest=1000;
		  ktest=-1000;
		  //cout<<"broken"<<endl;
		}
	      }
	    }
	  }
	}	     
      }
    }
	  
    cout<<" Second Tx Cor are: dX="<<DummyTx[0]-FinalTxCor[0]<<" ,dY="<<DummyTx[1]-FinalTxCor[1]<<" ,dZ="<<DummyTx[2]-FinalTxCor[2]<<"  Xi="<<DummyTx[0]<<" ,Yi="<<DummyTx[1]<<" ,Zi="<<DummyTx[2]<<"   Xi="<<FinalTxCor[0]<<" ,Yi="<<FinalTxCor[1]<<" ,Zi="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
	  
  }


  
}
