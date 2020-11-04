#include "Interferometer.cc"

void RunInterferometerMultTx(){
  DeclareAntennaConfig();
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx]; ////Channel Hit Time
  double ChDRTime[TotalAntennasRx]; ////Channel Hit Time
  
  double FinalTxCor[3];
  double JitterNumber=5;

  TH1D *hx=new TH1D("","",200,-500,500);
  TH1D *hy=new TH1D("","",200,-500,500);
  TH1D *hz=new TH1D("","",200,-500,500);

  TH1D *hxZoom=new TH1D("","",100,-50,50);
  TH1D *hyZoom=new TH1D("","",100,-50,50);
  TH1D *hzZoom=new TH1D("","",100,-50,50);
  
  double FinalMinValue;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;
  
  ofstream aout("RFTxPositions.txt");
  
  for(double i=-500; i<501;i=i+100){
    for(double j=-500; j<501;j=j+100){

      if((i!=0 && j==0) || (i==0 && j!=0) || (i!=0 && j!=0)){
	
	for(double k=-1; k>-1001;k=k-50){
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
		aout<<countRecoFailed<<" "<<DummyTx[0]<<" "<<DummyTx[1]<<" "<<DummyTx[2]<<endl;
		countRecoFailed++;
	      }
	    }
	    hx->Fill(i-FinalTxCor[0]);	
	    hy->Fill(j-FinalTxCor[1]);
	    hz->Fill(k-FinalTxCor[2]);

	    hxZoom->Fill(i-FinalTxCor[0]);	
	    hyZoom->Fill(j-FinalTxCor[1]);
	    hzZoom->Fill(k-FinalTxCor[2]);
	    
	    CheckReco=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
	    if(CheckReco<3){
	      countRecoSuccess++;
	    }
	 
	    countTotal++;
	  }////Check Station Trigger
	}////k loop
      }////check Tx is not at zero
    }////j loop
  }////i loop

  cout<<"SuccessReco: "<<countRecoSuccess<<", FailedReco: "<<countRecoFailed<<", TotalReco: "<<countTotal<<endl;

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

  hxZoom->SetTitle("#DeltaX; X_{true} - X_{reco} (m);No. of Events");
  hyZoom->SetTitle("#DeltaY; Y_{true} - Y_{reco} (m);No. of Events");
  hzZoom->SetTitle("#DeltaZ; Z_{true} - z_{reco} (m);No. of Events");
  TCanvas *c12=new TCanvas("c12","c12");
  c12->Divide(3);
  c12->cd(1);
  hxZoom->Draw();
  c12->cd(2);
  hyZoom->Draw();
  c12->cd(3);
  hzZoom->Draw();
 
  
}
