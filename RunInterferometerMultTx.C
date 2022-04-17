#include "Interferometer.cc"

void RunInterferometerMultTx(){
  DeclareAntennaConfig();

  double ExpectedTimeJitter=2;// in ns
  double ExpectedPositionUncertainty=5;//in m

  bool RefineRecoResults=false;
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];  
  
  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  double FinalMinValue=0;
  int FinalMinValueBin=0;
  int Iterations;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;
  
  double InitialTxCor_ThPhR[3]={0,0,0};
  double InitialTxCor_XYZ[3]={0,0,0};
  double FinalTxCor_ThPhR[3]={0,0,0};
  double FinalTxCor_XYZ[3]={0,0,0};
  double TrueTxCor_ThPhR[3]={0,0,0};
  double TrueTxCor_XYZ[3]={0,0,0};
  int IsItBelowStation;
  double ArrivalDirection[3];
 
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
	
	for(double k=-10; k>-501;k=k-100){
	  double DummyTxCor[3]={i,j,k};	  

	  for(int ixyz=0;ixyz<3;ixyz++){
	    double RandNum = 0;//(RandNumIni->Rndm(ixyz)*2-1)*ExpectedPositionUncertainty;
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

	  if(DummyTxCor[2]<0){
	    Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
	  }else{
	    Interferometer::GenerateChHitTimeAndCheckHits_Air(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
	  }
	  Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTime);
	  
	  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
	  if(CheckStationTrigger==true){
	    cout<<" X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
	    cout<<" Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
	    
	    double A=1.78;
	    double B=-0.43;
	    double C=0.0132;    
	    IceRayTracing::SetA(A);
	    IceRayTracing::SetB(B);
	    IceRayTracing::SetC(C);      
	    //Check if the station was hit from above or below
	    IsItBelowStation=Interferometer::IsItAboveOrBelow(ChHitTime,IgnoreCh) ;

	    // double X2a=Interferometer::GetChiSquaredThPhR(TrueTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, IsItBelowStation,0,50);    
	    // cout<<"chi sqaured values are "<<X2a<<endl;
      
	    Interferometer::GetRecieveAngle(ChHitTime, IgnoreCh, ChSNR, 100, ArrivalDirection);

	    //Set up variables for interferometry
	    double GuessResultCor[3][4]; 
    
	    //Start counting time it takes to guess the initial condition
	    auto t1 = std::chrono::high_resolution_clock::now();

	    ////Get a guess direction
	    Interferometer::GetApproximateMinThPhR(GuessResultCor,ChHitTime,IgnoreCh,ChSNR,IsItBelowStation,10);
      
	    //Get a guess distance
	    Interferometer::GetApproximateDistance(GuessResultCor,ChHitTime,IgnoreCh,ChSNR,IsItBelowStation,10);
      
	    //Calculate the time it took to guess the initial condition
	    auto t2 = std::chrono::high_resolution_clock::now();
	    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	    DurationInitialCondition=duration/1000;      
       
	    //Convert initial guess to radians
	    InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
	    InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
	    InitialTxCor_ThPhR[2]=GuessResultCor[0][2];  
      
	    //Set the radial width around guess distance value in meters    
	    double MinimizerRadialWidth=500;

	    //do the interferometry
	    Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   

	    if(RefineRecoResults==true){
	      cout<<"1st try Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]*(180./Interferometer::pi)<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]*(180./Interferometer::pi)<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
      
	      InitialTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
	      InitialTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180);
	      InitialTxCor_ThPhR[2]=FinalTxCor_ThPhR[2];  
	
	      Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   
	    }
	  
	    //Calcualte the total duration
	    DurationTotal=DurationInitialCondition+DurationReconstruction;

	    //Convert Initial guess direction from theta,phi,r to x,y,z
	    InitialTxCor_XYZ[0]=0;
	    InitialTxCor_XYZ[1]=0;
	    InitialTxCor_XYZ[2]=0;
	    Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
	    InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
	    InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 
      
	    //Convert final reconstructed direction from theta,phi,r to x,y,z
	    FinalTxCor_XYZ[0]=0;
	    FinalTxCor_XYZ[1]=0;
	    FinalTxCor_XYZ[2]=0;
	    FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
	    FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180); 
	    Interferometer::ThPhRtoXYZ(FinalTxCor_ThPhR, FinalTxCor_XYZ);
	    FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(180./Interferometer::pi);
	    FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(180./Interferometer::pi);
  	    
	    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
	    cout<<"True values are: |  Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
	    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<DurationTotal<<endl;	
	    cout<<" "<<endl;
	    
	    double dX=TrueTxCor_XYZ[0]-FinalTxCor_XYZ[0];
	    double dY=TrueTxCor_XYZ[1]-FinalTxCor_XYZ[1];
	    double dZ=TrueTxCor_XYZ[2]-FinalTxCor_XYZ[2];
	    
	    double CheckReco=sqrt(dX*dX+dY*dY+dZ*dZ);
	    if(CheckReco>ExpectedPositionUncertainty*2){
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

	    hDuration->Fill(DurationTotal);
	    hIterations->Fill(Iterations);
	    hDurationIterations->Fill(DurationTotal,Iterations);    
	 
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
