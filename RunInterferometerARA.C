#include "Interferometer.cc"

void GetCPCor(Int_t iStation, vector <double> CalPulCor[3],double year){

  ///////////////////////////////////////////////////////////////
  Double_t antLocD5[2][3];
  Double_t antLocD6[2][3];
  
  AraStationInfo *stationInfo=new AraStationInfo(iStation,year);
  AraGeomTool *geom = AraGeomTool::Instance();;
  
  Int_t numCalAnts=stationInfo->getNumCalAnts();
  for(int i=0;i<numCalAnts;i++) {
    AraCalAntennaInfo *calAntInfo=stationInfo->getCalAntennaInfo(i);
    //cout << calAntInfo->antName <<" "<<numCalAnts<<endl;
    //stationInfo->getCalAntennaInfo(i)->printAntennaInfo();
    Double_t *antLocTx=stationInfo->getCalAntennaInfo(i)->getLocationXYZ();
    /////////0 is H
    /////////1 is V
    if(i<2){
      antLocD5[i][0]=antLocTx[0];       
      antLocD5[i][1]=antLocTx[1];       
      antLocD5[i][2]=antLocTx[2];
      //cout<<i<<" D5: "<<antLocD5[i][0]<<" "<<antLocD5[i][1]<<" "<<antLocD5[i][2]<<endl;       
    }
    if(i>1){
      antLocD6[i-2][0]=antLocTx[0];       
      antLocD6[i-2][1]=antLocTx[1];       
      antLocD6[i-2][2]=antLocTx[2];       
      //cout<<i<<" D6: "<<antLocD6[i-2][0]<<" "<<antLocD6[i-2][1]<<" "<<antLocD6[i-2][2]<<endl;         
    }
    Double_t ThPhR[3]={0,0,0};
    Interferometer::XYZtoThPhR(antLocTx,ThPhR);
    ThPhR[0]=ThPhR[0]*(180./Interferometer::pi);
    ThPhR[1]=ThPhR[1]*(180./Interferometer::pi); 
    CalPulCor[0].push_back(ThPhR[0]);
    CalPulCor[1].push_back(ThPhR[1]);
    CalPulCor[2].push_back(ThPhR[2]);
    
  }
}

void RunInterferometerARA(){
  DeclareAntennaConfigARA(2);

  double ExpectedTimeJitter=10;// in ns
  double ExpectedPositionUncertainty=5;// in m
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];  
  
  double FinalMinValue=0;
  int FinalMinValueBin=0;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;

  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  vector <double> DurationReco;
  vector <double> DurationTot;
  vector <double> DurationIniCon;
 
  double InitialTxCor_ThPhR[3]={0,0,0};
  double InitialTxCor_XYZ[3]={0,0,0};
  double FinalTxCor_ThPhR[3]={0,0,0};
  double FinalTxCor_XYZ[3]={0,0,0};
  
  TH1D *hDuration=new TH1D("","",200,0,60000);
  TH1D *hIterations=new TH1D("","",50,0,50);

  TH2D *hDurationIterations=new TH2D("","",200,0,60000,50,0,50);
  
  TH1D *hx=new TH1D("","",200,-200,200);
  TH1D *hy=new TH1D("","",200,-200,200);
  TH1D *hz=new TH1D("","",200,-200,200);

  TH1D *hr=new TH1D("","",200,-100,100);
  TH1D *hth=new TH1D("","",200,-100,100);
  TH1D *hph=new TH1D("","",200,-100,100);
  
  ofstream aout("RFTxPositions.txt");

  TRandom3 *RandNumIni = new TRandom3(0);

  vector <double> CalPulCor[3];
  GetCPCor(2, CalPulCor,2017);
  
  double TrueTxCor_ThPhR[3]={CalPulCor[0][3],CalPulCor[1][3],CalPulCor[2][3]};
  TrueTxCor_ThPhR[0]=TrueTxCor_ThPhR[0]*(Interferometer::pi/180);
  TrueTxCor_ThPhR[1]=TrueTxCor_ThPhR[1]*(Interferometer::pi/180); 
  double TrueTxCor_XYZ[3]={0,0,0}; 
  
  Interferometer::ThPhRtoXYZ(TrueTxCor_ThPhR,TrueTxCor_XYZ);
  
  double DummyTxCor[3]={TrueTxCor_XYZ[0],TrueTxCor_XYZ[1],TrueTxCor_XYZ[2]};
  
  TrueTxCor_ThPhR[0]=TrueTxCor_ThPhR[0]*(180./Interferometer::pi);
  TrueTxCor_ThPhR[1]=TrueTxCor_ThPhR[1]*(180./Interferometer::pi);  

  double MinimizerRadialWidth=100;
  
  for(double k=0; k<400;k++){

    for(int ixyz=0;ixyz<3;ixyz++){
      double RandNum = (RandNumIni->Rndm(ixyz)*2-1)*ExpectedPositionUncertainty;
      InitialTxCor_XYZ[ixyz]=DummyTxCor[ixyz]+RandNum;
    }
    
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      for(int iray=0;iray<2;iray++){
	IgnoreCh[iray][iRx]=1;
	ChSNR[iray][iRx]=10;
	if(iRx>7 || iray==1){
	  IgnoreCh[iray][iRx]=0;
	}
      }
    }
    
    Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTime,IgnoreCh);
    Interferometer::AddGaussianJitterToHitTimes(ExpectedTimeJitter,ChHitTime);
    
    //bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
    //if(CheckStationTrigger==true){
   
      double GuessResultCor[3][3];
      Interferometer::GetApproximateMinThPhR(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
      Interferometer::GetApproximateDistance(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
      
      InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[2]=GuessResultCor[0][2];
      Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
      
      double Duration=0;
      int Iterations=0;  
      Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ExpectedPositionUncertainty, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, Duration, Iterations,MinimizerRadialWidth); 
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
      
      //cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<" | X_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]<<endl;
      cout<<k<<" Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<" | Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
      cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;
      
      double dX=TrueTxCor_XYZ[0]-FinalTxCor_XYZ[0];
      double dY=TrueTxCor_XYZ[1]-FinalTxCor_XYZ[1];
      double dZ=TrueTxCor_XYZ[2]-FinalTxCor_XYZ[2];

      double dTh=TrueTxCor_ThPhR[0]-FinalTxCor_ThPhR[0];
      double dPh=TrueTxCor_ThPhR[1]-FinalTxCor_ThPhR[1];
      double dR=TrueTxCor_ThPhR[2]-FinalTxCor_ThPhR[2];
      
      double CheckReco=sqrt(dX*dX+dY*dY+dZ*dZ);
      if(CheckReco>10){
	cout<<"Minimization Failed!"<<endl;
	cout<<countRecoFailed<<" "<<TrueTxCor_XYZ[0]<<" "<<TrueTxCor_XYZ[1]<<" "<<TrueTxCor_XYZ[2]<<" "<<CheckReco<<" "<<ExpectedPositionUncertainty<<endl;
	countRecoFailed++;
      }else{
	countRecoSuccess++;
      }
      hx->Fill(dX);	
      hy->Fill(dY);
      hz->Fill(dZ);
      
      hr->Fill(dR);	
      hth->Fill(dTh);
      hph->Fill(dPh);
      
      hDuration->Fill(Duration);
      hIterations->Fill(Iterations);
      hDurationIterations->Fill(Duration,Iterations);    
      
      countTotal++;
      //}////Check Station Trigger
  }////k loop
      
  cout<<"SuccessReco: "<<countRecoSuccess<<", FailedReco: "<<countRecoFailed<<", TotalReco: "<<countTotal<<endl;

  hx->SetTitle("#DeltaX; X_{true} - X_{reco} (m);No. of Events");
  hy->SetTitle("#DeltaY; Y_{true} - Y_{reco} (m);No. of Events");
  hz->SetTitle("#DeltaZ; Z_{true} - Z_{reco} (m);No. of Events");
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(3);
  c1->cd(1);
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  c1->cd(1)->SetLogy();
  hx->Draw();
  c1->cd(2);
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  c1->cd(2)->SetLogy();
  hy->Draw();
  c1->cd(3);
  c1->cd(3)->SetGridx();
  c1->cd(3)->SetGridy();
  c1->cd(3)->SetLogy();
  hz->Draw();

  hr->SetTitle("#DeltaR; R_{true} - R_{reco} (m);No. of Events");
  hth->SetTitle("#DeltaTheta; Theta_{true} - Theta_{reco} (m);No. of Events");
  hph->SetTitle("#DeltaPhi; Phi_{true} - Phi_{reco} (m);No. of Events");
  TCanvas *c12=new TCanvas("c12","c12");
  c12->Divide(3);
  c12->cd(1);
  c12->cd(1)->SetGridx();
  c12->cd(1)->SetGridy();
  c12->cd(1)->SetLogy();
  hth->Draw();
  c12->cd(2);
  c12->cd(2)->SetGridx();
  c12->cd(2)->SetGridy();
  c12->cd(2)->SetLogy();
  hph->Draw();
  c12->cd(3);
  c12->cd(3)->SetGridx();
  c12->cd(3)->SetGridy();
  c12->cd(3)->SetLogy();
  hr->Draw();
  
  TCanvas *c2=new TCanvas("c2","c2");
  c2->Divide(2);
  c2->cd(1);
  c2->cd(1)->SetGridx();
  c2->cd(1)->SetGridy();
  c2->cd(1)->SetLogy();
  hDuration->Draw();
  hDuration->SetTitle(";No. of Events;Reconstruction Duration (ms);");
  c2->cd(2);
  c2->cd(2)->SetGridx();
  c2->cd(2)->SetGridy();
  c2->cd(2)->SetLogy();
  hIterations->Draw();
  hIterations->SetTitle(";No. of Events;Total No. of Minimizer Iterations;");

  TCanvas *c3=new TCanvas("c3","c3");
  c3->cd(1);
  hDurationIterations->Draw("colz");
  hDurationIterations->SetTitle("Total No. of Minimizer Iterations;No. of Events");
  
}