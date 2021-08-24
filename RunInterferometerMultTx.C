#include "Interferometer.cc"

void RunInterferometerMultTx(){
  DeclareAntennaConfig();

  double ExpectedTimeJitter=2;// in ns
  double ExpectedPositionUncertainty=sqrt(5*5+pow(ExpectedTimeJitter*0.3*1.78,2));// in m
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];  
  
  double FinalMinValue=0;
  int FinalMinValueBin=0;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;
  
  double InitialTxCor_ThPhR[3]={0,0,0};
  double InitialTxCor_XYZ[3]={0,0,0};
  double FinalTxCor_ThPhR[3]={0,0,0};
  double FinalTxCor_XYZ[3]={0,0,0};
  double TrueTxCor_ThPhR[3]={0,0,0};
  double TrueTxCor_XYZ[3]={0,0,0};
  
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
	
	for(double k=-1; k>-501;k=k-100){
	  double DummyTxCor[3]={i,j,k};	  

	  if(DummyTxCor[2]>0){
	    DummyTxCor[2]=-1;
	  }
	  
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
	    cout<<"XYZ Tx coordinates are  "<<i<<" "<<j<<" "<<k<<endl;
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
  
	    cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<" | X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
	    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<" | Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
	    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;

	    double dX=TrueTxCor_XYZ[0]-FinalTxCor_XYZ[0];
	    double dY=TrueTxCor_XYZ[1]-FinalTxCor_XYZ[1];
	    double dZ=TrueTxCor_XYZ[2]-FinalTxCor_XYZ[2];
	    
	    double CheckReco=sqrt(dX*dX+dY*dY+dZ*dZ);
	    if(CheckReco>ExpectedPositionUncertainty){
	      cout<<"Minimization Failed!"<<endl;
	      cout<<countRecoFailed<<" "<<TrueTxCor_XYZ[0]<<" "<<TrueTxCor_XYZ[1]<<" "<<TrueTxCor_XYZ[2]<<" "<<CheckReco<<" "<<ExpectedPositionUncertainty<<endl;
	      countRecoFailed++;
	    }else{
	      countRecoSuccess++;
	    }
	    hx->Fill(i-FinalTxCor_XYZ[0]);	
	    hy->Fill(j-FinalTxCor_XYZ[1]);
	    hz->Fill(k-FinalTxCor_XYZ[2]);

	    hxZoom->Fill(i-FinalTxCor_XYZ[0]);	
	    hyZoom->Fill(j-FinalTxCor_XYZ[1]);
	    hzZoom->Fill(k-FinalTxCor_XYZ[2]);

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
