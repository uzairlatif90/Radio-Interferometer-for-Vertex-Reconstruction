#include "../../Interferometer/Interferometer.cc"

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

void PlotRecoResults(){
  
  Int_t Run=9129;

  int eventNum;
  double unixTime;
  double firstUnixTime;
  double timeStamp;
  int runNum;
  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  double FinalMinValue;
  int Iterations;
  bool isCalpulserTrig;
  bool isSoftwareTrig;
  double FinalTxCor_XYZ[3];
  double FinalTxCor_ThPhR[3];
  double FinalTxCor_XYZ_fR[3];
  double FinalTxCor_ThPhR_fR[3];
  double InitialTxCor_XYZ[3];
  double InitialTxCor_ThPhR[3];
  double VoltageSNR[16]; 
  double VoltageSNRV[3]; 
  int VoltageSNRV_Ch[3];
  double VoltageSNRH[3]; 
  int VoltageSNRH_Ch[3];
  
  double CorScore[16];
  double ChHitTime[2][16];
  int IgnoreCh[2][16];
 
  TString InputFileName="Run";
  InputFileName+=Run;
  //InputFileName+="AllEvents_varyingR_varyingIni.root";
  //InputFileName+="AllEvents.root";
  //InputFileName+="AllEvents_woCPini.root";
  InputFileName+="AllEvents_wCPini.root";  
  TFile *InputFile=new TFile(InputFileName);
  TString OutputFileName="Run";
  OutputFileName+=Run;
  OutputFileName+="RecoResults.root";
  TFile *OutputFile=new TFile(OutputFileName,"RECREATE");  
  
  TTree *RecoTree =(TTree*)InputFile->Get("RecoTree");

  RecoTree->SetBranchAddress("eventNum",&eventNum);
  RecoTree->SetBranchAddress("unixTime",&unixTime);
  RecoTree->SetBranchAddress("firstUnixTime",&firstUnixTime);
  RecoTree->SetBranchAddress("timeStamp",&timeStamp);  
  RecoTree->SetBranchAddress("runNum",&runNum);
  RecoTree->SetBranchAddress("DurationTotal",&DurationTotal);
  RecoTree->SetBranchAddress("DurationReconstruction",&DurationReconstruction);
  RecoTree->SetBranchAddress("DurationInitialCondition",&DurationInitialCondition);
  RecoTree->SetBranchAddress("FinalMinValue",&FinalMinValue);
  RecoTree->SetBranchAddress("Iterations",&Iterations);
  RecoTree->SetBranchAddress("isCalpulserTrig",&isCalpulserTrig);
  RecoTree->SetBranchAddress("isSoftwareTrig",&isSoftwareTrig);
  RecoTree->SetBranchAddress("FinalTxCor_XYZ",FinalTxCor_XYZ);
  RecoTree->SetBranchAddress("FinalTxCor_ThPhR",FinalTxCor_ThPhR);
  RecoTree->SetBranchAddress("FinalTxCor_XYZ_fR",FinalTxCor_XYZ_fR);
  RecoTree->SetBranchAddress("FinalTxCor_ThPhR_fR",FinalTxCor_ThPhR_fR);
  RecoTree->SetBranchAddress("InitialTxCor_XYZ",InitialTxCor_XYZ);
  RecoTree->SetBranchAddress("InitialTxCor_ThPhR",InitialTxCor_ThPhR);  
  RecoTree->SetBranchAddress("VoltageSNR",VoltageSNR); 
  // RecoTree->SetBranchAddress("VoltageSNRV",VoltageSNRV);
  // RecoTree->SetBranchAddress("VoltageSNRV_Ch",VoltageSNRV_Ch);
  // RecoTree->SetBranchAddress("VoltageSNRH",VoltageSNRH);
  // RecoTree->SetBranchAddress("VoltageSNRH_Ch",VoltageSNRH_Ch);
  RecoTree->SetBranchAddress("CorScore",CorScore);
  RecoTree->SetBranchAddress("ChHitTime",ChHitTime);
  RecoTree->SetBranchAddress("IgnoreCh",IgnoreCh); 
  
  double numEntries=RecoTree->GetEntries();

  RecoTree->GetEntry(0);

  TH1D *hIterations=new TH1D("","",200,0,50);
  
  TH1D * hReco_dXYZ[3];
  hReco_dXYZ[0]=new TH1D("","",200,-500,+500);
  hReco_dXYZ[1]=new TH1D("","",200,-500,+500);
  hReco_dXYZ[2]=new TH1D("","",200,-500,+500);

  TH1D * hReco_dThPhR[3];
  hReco_dThPhR[0]=new TH1D("","",200,-100,100);
  hReco_dThPhR[1]=new TH1D("","",200,-100,100);  
  hReco_dThPhR[2]=new TH1D("","",200,-100,100);

  TH1D * hReco_dThPhR_fR[3];
  hReco_dThPhR_fR[0]=new TH1D("","",200,-100,100);
  hReco_dThPhR_fR[1]=new TH1D("","",200,-100,100);  
  hReco_dThPhR_fR[2]=new TH1D("","",100,-50,50);  

  TH1D * hReco_ThPhR[3];
  hReco_ThPhR[0]=new TH1D("","",200,90,180);
  hReco_ThPhR[1]=new TH1D("","",200,-180,180);  
  hReco_ThPhR[2]=new TH1D("","",200,0,2000);

  TH1D * hReco_ThPhR_fR[3];
  hReco_ThPhR_fR[0]=new TH1D("","",200,90,180);
  hReco_ThPhR_fR[1]=new TH1D("","",200,-180,180);  
  hReco_ThPhR_fR[2]=new TH1D("","",200,0,2000);

  TH2D *hChHitTime=new TH2D("","",200,0,200,16,0,16);
  
  vector <double> IterationsReco;
  vector <double> DurationReco;
  vector <double> DurationTot;
  vector <double> DurationIniCon;
  vector <double> UnixTimeReco;
  vector <double> Reco_XYZ[3];
  vector <double> Reco_ThPhR[3];

  vector <double> CalPulCor[3];
  GetCPCor(2, CalPulCor,2017);

  int CPnum=0;
  if(Run==9129){
    CPnum=3;
  }
  if(Run==9129){
    CPnum=2;
  }
  
  double TrueTxCor_ThPhR[3]={CalPulCor[0][3],CalPulCor[1][3],CalPulCor[2][3]};
  TrueTxCor_ThPhR[0]=TrueTxCor_ThPhR[0]*(Interferometer::pi/180);
  TrueTxCor_ThPhR[1]=TrueTxCor_ThPhR[1]*(Interferometer::pi/180); 
  double TrueTxCor_XYZ[3]={0,0,0}; 
  
  Interferometer::ThPhRtoXYZ(TrueTxCor_ThPhR,TrueTxCor_XYZ);
  
  TrueTxCor_ThPhR[0]=TrueTxCor_ThPhR[0]*(180./Interferometer::pi);
  TrueTxCor_ThPhR[1]=TrueTxCor_ThPhR[1]*(180./Interferometer::pi);

  DeclareAntennaConfigARA(2);
  TLine * Lines[16];

  double ChHitTimeB[2][16];
  int IgnoreChB[2][16];
  for(Int_t ich=0; ich<16; ich++){
    IgnoreChB[0][ich]=1;
    if(Run==9129 && ich>7){
      IgnoreChB[0][ich]=0;
    }
    if(Run==9185 && ich<8){
      IgnoreChB[0][ich]=0;
    }
    IgnoreChB[1][ich]=0;
   
  }
 
  Interferometer::GenerateChHitTimeAndCheckHits(TrueTxCor_XYZ,ChHitTimeB,IgnoreChB);
  double ChDRTimeB[TotalAntennasRx];
  int ChHitOrderB[TotalAntennasRx];
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTimeB,IgnoreChB,ChDRTimeB,ChHitOrderB);

  for(Int_t ich=0; ich<16; ich++){
    if(IgnoreChB[0][ich]==1){
      Lines[ich]=new TLine(ChHitTimeB[0][ich],ich-0.5,ChHitTimeB[0][ich],ich+1.5);
      Lines[ich]->SetLineColor(kRed);
      Lines[ich]->SetLineWidth(3);
    }
  }
  
  ///Start looping over the events
  for(Int_t ievt=0; ievt<numEntries; ievt++){
 
    RecoTree->GetEntry(ievt);

    //if(InitialTxCor_XYZ[0]!=0 && InitialTxCor_XYZ[1]!=0 && InitialTxCor_XYZ[2]!=0){
    if(isCalpulserTrig==true){
      //if(eventNum==1301){
      //cout<<" its a calpulser event "<<eventNum<<endl;
      hIterations->Fill(Iterations);
      IterationsReco.push_back(Iterations);
      DurationReco.push_back(DurationReconstruction);
      DurationTot.push_back(DurationTotal);
      DurationIniCon.push_back(DurationInitialCondition);
      UnixTimeReco.push_back(unixTime-firstUnixTime);

      //cout<<setprecision(10)<<ievt<<" "<<unixTime<<" "<<firstUnixTime<<" "<<unixTime-firstUnixTime<<endl;
    
      Reco_XYZ[0].push_back(FinalTxCor_XYZ[0]);
      Reco_XYZ[1].push_back(FinalTxCor_XYZ[1]);
      Reco_XYZ[2].push_back(FinalTxCor_XYZ[2]);

      double dX=TrueTxCor_XYZ[0]-FinalTxCor_XYZ[0];
      double dY=TrueTxCor_XYZ[1]-FinalTxCor_XYZ[1];
      double dZ=TrueTxCor_XYZ[2]-FinalTxCor_XYZ[2];

      double dTh=TrueTxCor_ThPhR[0]-FinalTxCor_ThPhR[0];
      double dPh=TrueTxCor_ThPhR[1]-FinalTxCor_ThPhR[1];
      double dR=TrueTxCor_ThPhR[2]-FinalTxCor_ThPhR[2];

      double dTh_fR=TrueTxCor_ThPhR[0]-FinalTxCor_ThPhR_fR[0];
      double dPh_fR=TrueTxCor_ThPhR[1]-FinalTxCor_ThPhR_fR[1];
      double dR_fR=TrueTxCor_ThPhR[2]-FinalTxCor_ThPhR_fR[2];
      
      hReco_dXYZ[0]->Fill(dX);
      hReco_dXYZ[1]->Fill(dY);
      hReco_dXYZ[2]->Fill(dZ);

      hReco_dThPhR[0]->Fill(dTh);
      hReco_dThPhR[1]->Fill(dPh);
      hReco_dThPhR[2]->Fill(dR);

      hReco_dThPhR_fR[0]->Fill(dTh_fR);
      hReco_dThPhR_fR[1]->Fill(dPh_fR);
      hReco_dThPhR_fR[2]->Fill(dR_fR);
      
      Reco_ThPhR[0].push_back(FinalTxCor_ThPhR[0]);
      Reco_ThPhR[1].push_back(FinalTxCor_ThPhR[1]);
      Reco_ThPhR[2].push_back(FinalTxCor_ThPhR[2]);

      hReco_ThPhR[0]->Fill(FinalTxCor_ThPhR[0]);
      hReco_ThPhR[1]->Fill(FinalTxCor_ThPhR[1]);
      hReco_ThPhR[2]->Fill(FinalTxCor_ThPhR[2]);

      hReco_ThPhR_fR[0]->Fill(FinalTxCor_ThPhR_fR[0]);
      hReco_ThPhR_fR[1]->Fill(FinalTxCor_ThPhR_fR[1]);
      hReco_ThPhR_fR[2]->Fill(FinalTxCor_ThPhR_fR[2]);

      for(Int_t ich=0; ich<16; ich++){
	if(IgnoreCh[0][ich]!=0){
	  hChHitTime->Fill(ChHitTime[0][ich],ich);
	}
      }
      // hReco_ThPhR[0]->Fill(InitialTxCor_ThPhR[0]);
      // hReco_ThPhR[1]->Fill(InitialTxCor_ThPhR[1]);
      // hReco_ThPhR[2]->Fill(InitialTxCor_ThPhR[2]);
    
      // hReco_dThPhR[0]->Fill(dTh);
      // hReco_dThPhR[1]->Fill(dPh);
      // hReco_dThPhR[2]->Fill(dR);
    }
    
  }////event loop 

  TMultiGraph * mgReco_XYZ=new TMultiGraph();
  TGraph *grReco_XYZ[3];

  TMultiGraph * mgReco_ThPhR=new TMultiGraph();
  TGraph *grReco_ThPhR[3];
  
  for(Int_t ixyz=0; ixyz<3; ixyz++){
    grReco_XYZ[ixyz]=new TGraph(Reco_XYZ[ixyz].size(),UnixTimeReco.data(),Reco_XYZ[ixyz].data());
    grReco_XYZ[ixyz]->SetMarkerStyle(20);
    grReco_XYZ[ixyz]->SetMarkerColor(2+ixyz);

    grReco_ThPhR[ixyz]=new TGraph(Reco_ThPhR[ixyz].size(),UnixTimeReco.data(),Reco_ThPhR[ixyz].data());
    grReco_ThPhR[ixyz]->SetMarkerStyle(20);
    grReco_ThPhR[ixyz]->SetMarkerColor(2+ixyz);
  }

  TGraph *grDurationReco=new TGraph(DurationReco.size(),UnixTimeReco.data(),DurationReco.data());
  grDurationReco->SetMarkerStyle(20);
  grDurationReco->SetMarkerColor(kRed);
  grDurationReco->SetTitle("Reconstruction Time (ms)");  

  TGraph *grDurationTot=new TGraph(DurationTot.size(),UnixTimeReco.data(),DurationTot.data());
  grDurationTot->SetMarkerStyle(20);
  grDurationTot->SetMarkerColor(kBlue);
  grDurationTot->SetTitle("Total Time (ms)");
  
  TGraph *grDurationIniCon=new TGraph(DurationIniCon.size(),UnixTimeReco.data(),DurationIniCon.data());
  grDurationIniCon->SetMarkerStyle(20);
  grDurationIniCon->SetMarkerColor(kGreen);
  grDurationIniCon->SetTitle("Initial Condition Time (ms)");
  
  TGraph *grIterationsReco=new TGraph(IterationsReco.size(),UnixTimeReco.data(),IterationsReco.data());
  grIterationsReco->SetMarkerStyle(20);
  grIterationsReco->SetMarkerColor(2);

  TMultiGraph *mgDuration=new TMultiGraph();
  mgDuration->Add(grDurationTot);
  mgDuration->Add(grDurationIniCon);
  mgDuration->Add(grDurationReco);
  
  TCanvas *c1 = new TCanvas("c1","c1",1800,1800);
  c1->Divide(1,3);
  c1->cd(1);
  c1->cd(1)->SetGridx();
  c1->cd(1)->SetGridy();
  grReco_XYZ[0]->GetXaxis()->SetLabelSize(0.05);
  grReco_XYZ[0]->GetYaxis()->SetLabelSize(0.05);
  grReco_XYZ[0]->GetXaxis()->SetTitleSize(0.05);
  grReco_XYZ[0]->GetYaxis()->SetTitleSize(0.05);
  grReco_XYZ[0]->Draw("AP");
  grReco_XYZ[0]->SetTitle(";Unixtime (s);Reconstructed X (m);");
  OutputFile->cd();
  grReco_XYZ[0]->Write("grRecoX");
  
  c1->cd(2);
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  grReco_XYZ[1]->GetXaxis()->SetLabelSize(0.05);
  grReco_XYZ[1]->GetYaxis()->SetLabelSize(0.05);
  grReco_XYZ[1]->GetXaxis()->SetTitleSize(0.05);
  grReco_XYZ[1]->GetYaxis()->SetTitleSize(0.05);
  grReco_XYZ[1]->Draw("AP");
  grReco_XYZ[1]->SetTitle(";Unixtime (s);Reconstructed Y (m);");
  OutputFile->cd();
  grReco_XYZ[1]->Write("grRecoY");
  
  c1->cd(3);
  c1->cd(3)->SetGridx();
  c1->cd(3)->SetGridy();
  grReco_XYZ[2]->GetXaxis()->SetLabelSize(0.05);
  grReco_XYZ[2]->GetYaxis()->SetLabelSize(0.05);
  grReco_XYZ[2]->GetXaxis()->SetTitleSize(0.05);
  grReco_XYZ[2]->GetYaxis()->SetTitleSize(0.05);
  grReco_XYZ[2]->Draw("AP");
  grReco_XYZ[2]->SetTitle(";Unixtime (s);Reconstructed Z (m);");  
  OutputFile->cd();
  grReco_XYZ[2]->Write("grRecoZ");

  TString PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecoXYZ";
  PlotFileName+=".png";
  c1->SaveAs(PlotFileName);
  int stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c1->SaveAs(PlotFileName);
  
  TCanvas *c2 = new TCanvas("c2","c2",1800,1800);
  c2->Divide(1,3);
  c2->cd(1);
  c2->cd(1)->SetGridx();
  c2->cd(1)->SetGridy();
  grReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.1);
  grReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.1);
  grReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.1);
  grReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.1);

  grReco_ThPhR[0]->Draw("AP");
  grReco_ThPhR[0]->SetTitle(";Unixtime (s);Reconstructed Theta (^{o});");
  OutputFile->cd();
  grReco_ThPhR[0]->Write("grRecoTh");
  
  c2->cd(2);
  c2->cd(2)->SetGridx();
  c2->cd(2)->SetGridy();
  grReco_ThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  grReco_ThPhR[1]->GetYaxis()->SetLabelSize(0.05);
  grReco_ThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  grReco_ThPhR[1]->GetYaxis()->SetTitleSize(0.05);
  grReco_ThPhR[1]->Draw("AP");
  grReco_ThPhR[1]->SetTitle(";Unixtime (s);Reconstructed Phi (^{o});");
  OutputFile->cd();
  grReco_ThPhR[1]->Write("grRecoPh");
  
  c2->cd(3);
  c2->cd(3)->SetGridx();
  c2->cd(3)->SetGridy();
  grReco_ThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  grReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.06);
  grReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  grReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.06);
  grReco_ThPhR[2]->Draw("AP");
  grReco_ThPhR[2]->SetTitle(";Unixtime (s);Reconstructed Displacement (m);");  
  OutputFile->cd();
  grReco_ThPhR[2]->Write("grRecoR");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecoThPhR";
  PlotFileName+=".png";
  c2->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  // c2->SaveAs(PlotFileName);
  
  TCanvas *c3 = new TCanvas("c3","c3",1800,1800);
  c3->Divide(2,3);
  c3->cd(1);
  c3->cd(1)->SetGridx();
  c3->cd(1)->SetGridy();
  c3->cd(1)->SetLogy();
  hReco_dXYZ[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dXYZ[0]->GetYaxis()->SetLabelSize(0.06);
  hReco_dXYZ[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dXYZ[0]->GetYaxis()->SetTitleSize(0.06);
  hReco_dXYZ[0]->Draw();
  hReco_dXYZ[0]->SetTitle(";Reconstructed X True - Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_dXYZ[0]->Write("hReco_dX");
  
  c3->cd(3);
  c3->cd(3)->SetGridx();
  c3->cd(3)->SetGridy();
  c3->cd(3)->SetLogy();
  hReco_dXYZ[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dXYZ[1]->GetYaxis()->SetLabelSize(0.06);
  hReco_dXYZ[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dXYZ[1]->GetYaxis()->SetTitleSize(0.06);
  hReco_dXYZ[1]->Draw();
  hReco_dXYZ[1]->SetTitle(";Reconstructed Y True - Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_dXYZ[1]->Write("hReco_dY");

  c3->cd(5);
  c3->cd(5)->SetGridx();
  c3->cd(5)->SetGridy();
  c3->cd(5)->SetLogy();
  hReco_dXYZ[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dXYZ[2]->GetYaxis()->SetLabelSize(0.06);
  hReco_dXYZ[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dXYZ[2]->GetYaxis()->SetTitleSize(0.06);
  hReco_dXYZ[2]->Draw();
  hReco_dXYZ[2]->SetTitle(";Reconstructed Z True - Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_dXYZ[2]->Write("hReco_dZ");

  c3->cd(2);
  c3->cd(2)->SetGridx();
  c3->cd(2)->SetGridy();
  //c3->cd(2)->SetLogy();
  hReco_dXYZ[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dXYZ[0]->GetYaxis()->SetLabelSize(0.06);
  hReco_dXYZ[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dXYZ[0]->GetYaxis()->SetTitleSize(0.06);
  hReco_dXYZ[0]->Draw();
  hReco_dXYZ[0]->SetTitle(";Reconstructed X True - Reco (m);No. of Events;");
  
  c3->cd(4);
  c3->cd(4)->SetGridx();
  c3->cd(4)->SetGridy();
  //c3->cd(4)->SetLogy();
  hReco_dXYZ[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dXYZ[1]->GetYaxis()->SetLabelSize(0.06);
  hReco_dXYZ[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dXYZ[1]->GetYaxis()->SetTitleSize(0.06);
  hReco_dXYZ[1]->Draw();
  hReco_dXYZ[1]->SetTitle(";Reconstructed Y True - Reco (m);No. of Events;");
 
  c3->cd(6);
  c3->cd(6)->SetGridx();
  c3->cd(6)->SetGridy();
  //c3->cd(6)->SetLogy();
  hReco_dXYZ[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dXYZ[2]->GetYaxis()->SetLabelSize(0.06);
  hReco_dXYZ[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dXYZ[2]->GetYaxis()->SetTitleSize(0.06);
  hReco_dXYZ[2]->Draw();
  hReco_dXYZ[2]->SetTitle(";Reconstructed Z True - Reco (m);No. of Events;");
  
  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecodXYZ";
  PlotFileName+=".png";
  c3->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c3->SaveAs(PlotFileName);

  TCanvas *c4 = new TCanvas("c4","c4",1800,1800);
  c4->Divide(2,3);
  c4->cd(1);
  c4->cd(1)->SetGridx();
  c4->cd(1)->SetGridy();
  c4->cd(1)->SetLogy();
  hReco_dThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetLabelSize(0.05);
  hReco_dThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetTitleSize(0.05);
  hReco_dThPhR[0]->Draw();
  hReco_dThPhR[0]->SetTitle(";Reconstructed Theta True - Reco (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR[0]->Write("hReco_dTh");
  
  c4->cd(3);
  c4->cd(3)->SetGridx();
  c4->cd(3)->SetGridy();
  c4->cd(3)->SetLogy();
  hReco_dThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[1]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[1]->GetYaxis()->SetTitleSize(0.06);
  hReco_dThPhR[1]->Draw();
  hReco_dThPhR[1]->SetTitle(";Reconstructed Phi True - Reco (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR[1]->Write("hReco_dPh");

  c4->cd(5);
  c4->cd(5)->SetGridx();
  c4->cd(5)->SetGridy();
  c4->cd(5)->SetLogy();
  hReco_dThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[2]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[2]->GetYaxis()->SetTitleSize(0.06);
  hReco_dThPhR[2]->Draw();
  hReco_dThPhR[2]->SetTitle(";Reconstructed Displacement True - Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_dThPhR[2]->Write("hReco_dR");

  c4->cd(2);
  c4->cd(2)->SetGridx();
  c4->cd(2)->SetGridy();
  hReco_dThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetTitleSize(0.06);
  hReco_dThPhR[0]->Draw();
  hReco_dThPhR[0]->SetTitle(";Reconstructed Theta True - Reco (^{o});No. of Events;");
  
  c4->cd(4);
  c4->cd(4)->SetGridx();
  c4->cd(4)->SetGridy();
  hReco_dThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[1]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[1]->GetYaxis()->SetTitleSize(0.06);
  hReco_dThPhR[1]->Draw();
  hReco_dThPhR[1]->SetTitle(";Reconstructed Phi True - Reco (^{o});No. of Events;");
 
  c4->cd(6);
  c4->cd(6)->SetGridx();
  c4->cd(6)->SetGridy();
  hReco_dThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[2]->GetYaxis()->SetLabelSize(0.05);
  hReco_dThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[2]->GetYaxis()->SetTitleSize(0.05);
  hReco_dThPhR[2]->Draw();
  hReco_dThPhR[2]->SetTitle(";Reconstructed Displacement True - Reco (m);No. of Events;");
 
  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecodThPhR";
  PlotFileName+=".png";
  c4->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c4->SaveAs(PlotFileName);

  TCanvas *c4b = new TCanvas("c4b","c4b",1800,1800);
  c4b->Divide(1,3);
  c4b->cd(1);
  c4b->cd(1)->SetGridx();
  c4b->cd(1)->SetGridy();
  c4b->cd(1)->SetLogy();
  hReco_dThPhR_fR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_fR[0]->GetYaxis()->SetLabelSize(0.05);
  hReco_dThPhR_fR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_fR[0]->GetYaxis()->SetTitleSize(0.05);
  hReco_dThPhR_fR[0]->Draw();
  hReco_dThPhR_fR[0]->SetTitle(";Reconstructed Theta True - Reco (Fixed R) (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_fR[0]->Write("hReco_dTh_fR");
  
  c4b->cd(2);
  c4b->cd(2)->SetGridx();
  c4b->cd(2)->SetGridy();
  c4b->cd(2)->SetLogy();
  hReco_dThPhR_fR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_fR[1]->GetYaxis()->SetLabelSize(0.05);
  hReco_dThPhR_fR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_fR[1]->GetYaxis()->SetTitleSize(0.05);
  hReco_dThPhR_fR[1]->Draw();
  hReco_dThPhR_fR[1]->SetTitle(";Reconstructed Phi True - Reco (Fixed R)(^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_fR[1]->Write("hReco_dPh_fR");

  c4b->cd(3);
  c4b->cd(3)->SetGridx();
  c4b->cd(3)->SetGridy();
  c4b->cd(3)->SetLogy();
  hReco_dThPhR_fR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_fR[2]->GetYaxis()->SetLabelSize(0.05);
  hReco_dThPhR_fR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_fR[2]->GetYaxis()->SetTitleSize(0.05);
  hReco_dThPhR_fR[2]->Draw();
  hReco_dThPhR_fR[2]->SetTitle(";Reconstructed Displacement True - Reco (Fixed R)(m);No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_fR[2]->Write("hReco_dR_fR");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecodThPhR_fR";
  PlotFileName+=".png";
  c4b->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c4b->SaveAs(PlotFileName);
  
  TCanvas *c5 = new TCanvas("c5","c5",1800,1800);
  c5->cd(1);
  c5->cd(1)->SetGridx();
  c5->cd(1)->SetGridy();
  mgDuration->GetXaxis()->SetLabelSize(0.04);
  mgDuration->GetYaxis()->SetLabelSize(0.04);
  mgDuration->GetXaxis()->SetTitleSize(0.04);
  mgDuration->GetYaxis()->SetTitleSize(0.04);
  mgDuration->Draw("AP");
  mgDuration->GetXaxis()->SetTitle("Run Unixtime (s)");
  mgDuration->GetYaxis()->SetTitle("Time Taken (ms)");
  c5->cd(1)->BuildLegend();
  OutputFile->cd();
  mgDuration->Write("mgDuration");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="Durations";
  PlotFileName+=".png";
  c5->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c5->SaveAs(PlotFileName);

  TCanvas *c6 = new TCanvas("c6","c6",1800,1800);
  c6->cd(1);
  c6->cd(1)->SetGridx();
  c6->cd(1)->SetGridy();
  grIterationsReco->GetXaxis()->SetLabelSize(0.05);
  grIterationsReco->GetYaxis()->SetLabelSize(0.05);
  grIterationsReco->GetXaxis()->SetTitleSize(0.05);
  grIterationsReco->GetYaxis()->SetTitleSize(0.05);
  grIterationsReco->Draw("AP");
  grIterationsReco->SetTitle(";Run Unixtime (s);Iterations taken by Minimizer to do Reco.;");  
  OutputFile->cd();
  grIterationsReco->Write("grRecoIterations");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecoIterations";
  PlotFileName+=".png";
  c6->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c6->SaveAs(PlotFileName);

  TCanvas *c7 = new TCanvas("c7","c7",1800,1800);
  c7->cd(1);
  c7->cd(1)->SetGridx();
  c7->cd(1)->SetGridy();
  hChHitTime->SetTitle("Channel Hit Time Pattern;Channel Hit Times (ns);Channel number; No. of Events;");
  hChHitTime->Draw("colz");
  for(Int_t ich=0; ich<16; ich++){
    if(IgnoreChB[0][ich]==1){
      Lines[ich]->Draw();
    }
  }
  OutputFile->cd();
  hChHitTime->Write("hChHitTime");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="hChHitTime";
  PlotFileName+=".png";
  c7->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c7->SaveAs(PlotFileName);
 
  TCanvas *c8 = new TCanvas("c8","c8",1800,1800);
  c8->Divide(1,3);
  c8->cd(1);
  c8->cd(1)->SetGridx();
  c8->cd(1)->SetGridy();
  c8->cd(1)->SetLogy();
  hReco_ThPhR[0]->GetXaxis()->SetNdivisions(20);
  hReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.05);
  hReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.05);
  hReco_ThPhR[0]->Draw();
  hReco_ThPhR[0]->SetTitle(";Reconstructed Theta Reco (^{o});No. of Events;");
  OutputFile->cd();
  hReco_ThPhR[0]->Write("hReco_Th");
  
  c8->cd(2);
  c8->cd(2)->SetGridx();
  c8->cd(2)->SetGridy();
  c8->cd(2)->SetLogy();
  hReco_ThPhR[1]->GetXaxis()->SetNdivisions(20);
  hReco_ThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR[1]->GetYaxis()->SetLabelSize(0.05);
  hReco_ThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR[1]->GetYaxis()->SetTitleSize(0.05);
  hReco_ThPhR[1]->Draw();
  hReco_ThPhR[1]->SetTitle(";Reconstructed Phi Reco (^{o});No. of Events;");
  OutputFile->cd();
  hReco_ThPhR[1]->Write("hReco_Ph");

  c8->cd(3);
  c8->cd(3)->SetGridx();
  c8->cd(3)->SetGridy();
  c8->cd(3)->SetLogy();
  hReco_ThPhR[2]->GetXaxis()->SetNdivisions(20);
  hReco_ThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.05);
  hReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.05);
  hReco_ThPhR[2]->Draw();
  hReco_ThPhR[2]->SetTitle(";Reconstructed Displacement Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_ThPhR[2]->Write("hReco_R");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="hRecoThPhR";
  PlotFileName+=".png";
  c8->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c8->SaveAs(PlotFileName);

  TCanvas *c8b = new TCanvas("c8b","c8b",1800,1800);
  c8b->Divide(1,3);
  c8b->cd(1);
  c8b->cd(1)->SetGridx();
  c8b->cd(1)->SetGridy();
  c8b->cd(1)->SetLogy();
  hReco_ThPhR_fR[0]->GetXaxis()->SetNdivisions(20);
  hReco_ThPhR_fR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR_fR[0]->GetYaxis()->SetLabelSize(0.05);
  hReco_ThPhR_fR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR_fR[0]->GetYaxis()->SetTitleSize(0.05);
  hReco_ThPhR_fR[0]->Draw();
  hReco_ThPhR_fR[0]->SetTitle(";Reconstructed Theta Reco (Fixed R) (^{o});No. of Events;");
  OutputFile->cd();
  hReco_ThPhR_fR[0]->Write("hReco_Th_fR");
  
  c8b->cd(2);
  c8b->cd(2)->SetGridx();
  c8b->cd(2)->SetGridy();
  c8b->cd(2)->SetLogy();
  hReco_ThPhR_fR[1]->GetXaxis()->SetNdivisions(20);
  hReco_ThPhR_fR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR_fR[1]->GetYaxis()->SetLabelSize(0.05);
  hReco_ThPhR_fR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR_fR[1]->GetYaxis()->SetTitleSize(0.05);
  hReco_ThPhR_fR[1]->Draw();
  hReco_ThPhR_fR[1]->SetTitle(";Reconstructed Phi Reco (Fixed R) (^{o});No. of Events;");
  OutputFile->cd();
  hReco_ThPhR_fR[1]->Write("hReco_Ph_fR");

  c8b->cd(3);
  c8b->cd(3)->SetGridx();
  c8b->cd(3)->SetGridy();
  c8b->cd(3)->SetLogy();
  hReco_ThPhR_fR[2]->GetXaxis()->SetNdivisions(20);
  hReco_ThPhR_fR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR_fR[2]->GetYaxis()->SetLabelSize(0.05);
  hReco_ThPhR_fR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR_fR[2]->GetYaxis()->SetTitleSize(0.05);
  hReco_ThPhR_fR[2]->Draw();
  hReco_ThPhR_fR[2]->SetTitle(";Displacement Reco (Fixed R) (m);No. of Events;");
  OutputFile->cd();
  hReco_ThPhR_fR[2]->Write("hReco_R_fR");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="hRecoThPhR_fR";
  PlotFileName+=".png";
  c8b->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c8b->SaveAs(PlotFileName);
  
}
