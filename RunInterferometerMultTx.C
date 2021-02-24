#include "Interferometer.cc"

void RunInterferometerMultTx(){
  DeclareAntennaConfig();

  double ExpectedTimeJitter=5;// in ns
  double ExpectedUncertaintyPosition=sqrt(5*5+pow(ExpectedTimeJitter*0.3*1.78,2));// in m
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];  
  
  double FinalMinValue=0;
  int FinalMinValueBin=0;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;
  
  double InitialTxCor[3];
  double ExpectedTxCor[3];
  double FinalTxCor[3];

  TH1D *hDuration=new TH1D("","",100,0,4000);
  TH1D *hIterations=new TH1D("","",100,0,500);

  TH2D *hDurationIterations=new TH2D("","",100,0,4000,100,0,500);
  
  TH1D *hx=new TH1D("","",200,-500,500);
  TH1D *hy=new TH1D("","",200,-500,500);
  TH1D *hz=new TH1D("","",200,-500,500);

  TH1D *hxZoom=new TH1D("","",100,-50,50);
  TH1D *hyZoom=new TH1D("","",100,-50,50);
  TH1D *hzZoom=new TH1D("","",100,-50,50);
  
  ofstream aout("RFTxPositions.txt");

  TRandom3 *RandNumIni = new TRandom3(0);
  
  for(double i=-500; i<501;i=i+100){
    for(double j=-500; j<501;j=j+100){

      if((i!=0 && j==0) || (i==0 && j!=0) || (i!=0 && j!=0)){
	
	for(double k=-1; k>-1001;k=k-50){
	  double DummyTxCor[3]={i,j,k};	  

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
	    CheckReco=sqrt(pow(i-FinalTxCor[0],2)+pow(j-FinalTxCor[1],2)+pow(k-FinalTxCor[2],2))/3;
	    cout<<"Final Reco Results are: dX="<<dX<<" ,dY="<<dY<<" ,dZ="<<dZ<<" |  Xtrue="<<DummyTxCor[0]<<" ,Ytrue="<<DummyTxCor[1]<<" ,Ztrue="<<DummyTxCor[2]<<" | Xreco="<<FinalTxCor[0]<<" ,Yreco="<<FinalTxCor[1]<<" ,Zreco="<<FinalTxCor[2]<<" | Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;
	    if(CheckReco>ExpectedUncertaintyPosition){
	      cout<<"Minimization Failed!"<<endl;
	      aout<<countRecoFailed<<" "<<ExpectedTxCor[0]<<" "<<ExpectedTxCor[1]<<" "<<ExpectedTxCor[2]<<endl;
	      countRecoFailed++;
	    }else{
	      countRecoSuccess++;
	    }
	    hx->Fill(i-FinalTxCor[0]);	
	    hy->Fill(j-FinalTxCor[1]);
	    hz->Fill(k-FinalTxCor[2]);

	    hxZoom->Fill(i-FinalTxCor[0]);	
	    hyZoom->Fill(j-FinalTxCor[1]);
	    hzZoom->Fill(k-FinalTxCor[2]);

	    hDuration->Fill(Duration);
	    hIterations->Fill(Iterations);
	    hDurationIterations->Fill(Duration,Iterations);    
	 
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

  TCanvas *c2=new TCanvas("c2","c2");
  c2->Divide(2);
  c2->cd(1);
  hDuration->Draw();
  hDuration->SetTitle("Reconstruction Duration (ms);No. of Events");
  c2->cd(2);
  hIterations->Draw();
  hIterations->SetTitle("Total No. of Minimizer Iterations;No. of Events");

  TCanvas *c3=new TCanvas("c3","c3");
  c3->cd(1);
  hDurationIterations->Draw("colz");
  hDurationIterations->SetTitle("Total No. of Minimizer Iterations;No. of Events");
  
}
