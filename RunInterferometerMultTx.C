#include "Interferometer.cc"

void RunInterferometerMultTx(){
  DeclareAntennaConfig();
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx]; ////Channel Hit Time

  double FinalTxCor[3];
  double JitterNumber=2;

  TH1D *hx=new TH1D("","",200,-500,500);
  TH1D *hy=new TH1D("","",200,-500,500);
  TH1D *hz=new TH1D("","",200,-500,500);

  double FinalMinValue;
  int count=0;
  int count2=0;
  int count3=0;
  int count4=0;
  int count5=0;
  
  for(double i=-500; i<501;i=i+100){
    for(double j=-500; j<501;j=j+100){

      if(i!=0 || j!=0){
	
	for(double k=-1; k>-1001;k=k-50){
	  double DummyTx[3]={i,j,k};

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
	    //cout<<"First Tx Cor are: dX="<<i-FinalTxCor[0]<<" ,dY="<<j-FinalTxCor[1]<<" ,dZ="<<k-FinalTxCor[2]<<"  Xi="<<i<<" ,Yi="<<j<<" ,Zi="<<k<<"   Xi="<<FinalTxCor[0]<<" ,Yi="<<FinalTxCor[1]<<" ,Zi="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
	    double checkdiff=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
	    if(checkdiff>3){
	      Interferometer::FindRoots3D(DummyTx,FinalTxCor,ChHitTime,IgnoreCh,0,status);
	      count5++;
	      
	      if(status!=0){
		count2++;
		Interferometer::FindRoots3D(FinalTxCor,FinalTxCor,ChHitTime,IgnoreCh,1,status);
	      }
	    
	      if(status!=0){
		count4++;
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
	    hx->Fill(i-FinalTxCor[0]);	
	    hy->Fill(j-FinalTxCor[1]);
	    hz->Fill(k-FinalTxCor[2]);
	  
	    //cout<<" Second Tx Cor are: dX="<<i-FinalTxCor[0]<<" ,dY="<<j-FinalTxCor[1]<<" ,dZ="<<k-FinalTxCor[2]<<"  Xi="<<i<<" ,Yi="<<j<<" ,Zi="<<k<<"   Xi="<<FinalTxCor[0]<<" ,Yi="<<FinalTxCor[1]<<" ,Zi="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
	    checkdiff=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
	    if(checkdiff<3){
	      count++;
	    }

	    if(checkdiff>3 && status==0){
	      count3++;
	    }
	  
	  }
	}
      }
    }
  }

  cout<<"count is "<<count<<" "<<count2<<" "<<count3<<" "<<count4<<endl;

  hx->SetTitle("#DeltaX; X_{true} - X_{reco} (m);No. of Events");
  hy->SetTitle("#DeltaY; Y_{true} - Y_{reco} (m);No. of Events");
  hz->SetTitle("#DeltaZ; Z_{true} - z_{reco} (m);No. of Events");
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(3);
  c1->cd(1);
  hx->Draw();
  c1->cd(2);
  hy->Draw();
  c1->cd(3);
  hz->Draw();
  
}
