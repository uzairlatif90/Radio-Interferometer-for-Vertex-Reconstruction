void PlotRecoResults(){

  Int_t Run=12576;

  double eventNum;
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
  double dXYZ[3];
  double dThPhR[3];
  
  TString InputFileName="Run";
  InputFileName+=Run;
  InputFileName+="AllEvents.root";  
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
  RecoTree->SetBranchAddress("dXYZ",dXYZ);
  RecoTree->SetBranchAddress("dThPhR",dThPhR);

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

  TH1D * hReco_ThPhR_fR[3];
  hReco_ThPhR_fR[0]=new TH1D("","",200,90,180);
  hReco_ThPhR_fR[1]=new TH1D("","",200,-180,180);  
  hReco_ThPhR_fR[2]=new TH1D("","",200,0,2000);
  
  vector <double> IterationsReco;
  vector <double> DurationReco;
  vector <double> DurationTot;
  vector <double> DurationIniCon;
  vector <double> UnixTimeReco;
  vector <double> Reco_XYZ[3];
  vector <double> Reco_ThPhR[3];

  vector <double> True_XYZ[3];
  vector <double> True_ThPhR[3];

  
  ///Start looping over the events
  for(Int_t ievt=0; ievt<numEntries; ievt++){
 
    RecoTree->GetEntry(ievt);
    if(InitialTxCor_XYZ[0]!=0){
      hIterations->Fill(Iterations);
      IterationsReco.push_back(Iterations);
      DurationReco.push_back(DurationReconstruction/1000);
      DurationTot.push_back(DurationTotal/1000);
      DurationIniCon.push_back(DurationInitialCondition/1000);
      UnixTimeReco.push_back(unixTime-firstUnixTime);

      //cout<<setprecision(10)<<ievt<<" "<<unixTime<<" "<<firstUnixTime<<" "<<unixTime-firstUnixTime<<endl;
    
      // Reco_XYZ[0].push_back(FinalTxCor_XYZ_fR[0]);
      // Reco_XYZ[1].push_back(FinalTxCor_XYZ_fR[1]);
      // Reco_XYZ[2].push_back(FinalTxCor_XYZ_fR[2]);
      // hReco_dXYZ[0]->Fill(InitialTxCor_XYZ[0]-FinalTxCor_XYZ_fR[0]);
      // hReco_dXYZ[1]->Fill(InitialTxCor_XYZ[1]-FinalTxCor_XYZ_fR[1]);
      // hReco_dXYZ[2]->Fill(InitialTxCor_XYZ[2]-FinalTxCor_XYZ_fR[2]);

      // Reco_ThPhR[0].push_back(FinalTxCor_ThPhR_fR[0]);
      // Reco_ThPhR[1].push_back(FinalTxCor_ThPhR_fR[1]);
      // Reco_ThPhR[2].push_back(FinalTxCor_ThPhR_fR[2]);
      // hReco_dThPhR[0]->Fill(InitialTxCor_ThPhR[0]-FinalTxCor_ThPhR_fR[0]);
      // hReco_dThPhR[1]->Fill(InitialTxCor_ThPhR[1]-FinalTxCor_ThPhR_fR[1]);
      // hReco_dThPhR[2]->Fill(InitialTxCor_ThPhR[2]-FinalTxCor_ThPhR_fR[2]);

      Reco_XYZ[0].push_back(FinalTxCor_XYZ[0]);
      Reco_XYZ[1].push_back(FinalTxCor_XYZ[1]);
      Reco_XYZ[2].push_back(FinalTxCor_XYZ[2]);
      hReco_dXYZ[0]->Fill(InitialTxCor_XYZ[0]-FinalTxCor_XYZ[0]);
      hReco_dXYZ[1]->Fill(InitialTxCor_XYZ[1]-FinalTxCor_XYZ[1]);
      hReco_dXYZ[2]->Fill(InitialTxCor_XYZ[2]-FinalTxCor_XYZ[2]);

      Reco_ThPhR[0].push_back(FinalTxCor_ThPhR[0]);
      Reco_ThPhR[1].push_back(FinalTxCor_ThPhR[1]);
      Reco_ThPhR[2].push_back(FinalTxCor_ThPhR[2]);
      hReco_dThPhR[0]->Fill(InitialTxCor_ThPhR[0]-FinalTxCor_ThPhR[0]);
      hReco_dThPhR[1]->Fill(InitialTxCor_ThPhR[1]-FinalTxCor_ThPhR[1]);
      hReco_dThPhR[2]->Fill(InitialTxCor_ThPhR[2]-FinalTxCor_ThPhR[2]);

      
      True_XYZ[0].push_back(InitialTxCor_XYZ[0]);
      True_XYZ[1].push_back(InitialTxCor_XYZ[1]);
      True_XYZ[2].push_back(InitialTxCor_XYZ[2]);

      True_ThPhR[0].push_back(InitialTxCor_ThPhR[0]);
      True_ThPhR[1].push_back(InitialTxCor_ThPhR[1]);
      True_ThPhR[2].push_back(InitialTxCor_ThPhR[2]);
    }
  }////event loop

  TMultiGraph * mgReco_XYZ[3];
  mgReco_XYZ[0]=new TMultiGraph();
  mgReco_XYZ[1]=new TMultiGraph();
  mgReco_XYZ[2]=new TMultiGraph();
  TGraph *grReco_XYZ[3];
  TGraph *grTrue_XYZ[3];

  TMultiGraph * mgReco_ThPhR[3];
  mgReco_ThPhR[0]=new TMultiGraph();
  mgReco_ThPhR[1]=new TMultiGraph();
  mgReco_ThPhR[2]=new TMultiGraph();
  TGraph *grReco_ThPhR[3];
  TGraph *grTrue_ThPhR[3];
  
  for(Int_t ixyz=0; ixyz<3; ixyz++){
    grReco_XYZ[ixyz]=new TGraph(Reco_XYZ[ixyz].size(),UnixTimeReco.data(),Reco_XYZ[ixyz].data());
    grReco_XYZ[ixyz]->SetMarkerStyle(20);
    grReco_XYZ[ixyz]->SetMarkerColor(2+ixyz);

    grReco_ThPhR[ixyz]=new TGraph(Reco_ThPhR[ixyz].size(),UnixTimeReco.data(),Reco_ThPhR[ixyz].data());
    grReco_ThPhR[ixyz]->SetMarkerStyle(20);
    grReco_ThPhR[ixyz]->SetMarkerColor(2+ixyz);

    grTrue_XYZ[ixyz]=new TGraph(True_XYZ[ixyz].size(),UnixTimeReco.data(),True_XYZ[ixyz].data());
    grTrue_XYZ[ixyz]->SetMarkerStyle(43);
    grTrue_XYZ[ixyz]->SetMarkerColor(1);

    grTrue_ThPhR[ixyz]=new TGraph(True_ThPhR[ixyz].size(),UnixTimeReco.data(),True_ThPhR[ixyz].data());
    grTrue_ThPhR[ixyz]->SetMarkerStyle(43);
    grTrue_ThPhR[ixyz]->SetMarkerColor(1);

    mgReco_XYZ[ixyz]->Add(grReco_XYZ[ixyz]);
    mgReco_XYZ[ixyz]->Add(grTrue_XYZ[ixyz]);

    mgReco_ThPhR[ixyz]->Add(grReco_ThPhR[ixyz]);
    mgReco_ThPhR[ixyz]->Add(grTrue_ThPhR[ixyz]);
    
  }

  TGraph *grDurationReco=new TGraph(DurationReco.size(),UnixTimeReco.data(),DurationReco.data());
  grDurationReco->SetMarkerStyle(20);
  grDurationReco->SetMarkerColor(kRed);
  grDurationReco->SetTitle("Reconstruction Time (s)");  

  TGraph *grDurationTot=new TGraph(DurationTot.size(),UnixTimeReco.data(),DurationTot.data());
  grDurationTot->SetMarkerStyle(20);
  grDurationTot->SetMarkerColor(kBlue);
  grDurationTot->SetTitle("Total Time (s)");
  
  TGraph *grDurationIniCon=new TGraph(DurationIniCon.size(),UnixTimeReco.data(),DurationIniCon.data());
  grDurationIniCon->SetMarkerStyle(20);
  grDurationIniCon->SetMarkerColor(kGreen);
  grDurationIniCon->SetTitle("Initial Condition Time (s)");
  
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
  mgReco_XYZ[0]->GetXaxis()->SetLabelSize(0.05);
  mgReco_XYZ[0]->GetYaxis()->SetLabelSize(0.06);
  mgReco_XYZ[0]->GetXaxis()->SetTitleSize(0.05);
  mgReco_XYZ[0]->GetYaxis()->SetTitleSize(0.06);
  mgReco_XYZ[0]->GetYaxis()->SetTitleOffset(0.7);
  mgReco_XYZ[0]->Draw("AP");
  mgReco_XYZ[0]->SetTitle(";Unixtime (s);Reconstructed X (m);");
  OutputFile->cd();
  mgReco_XYZ[0]->Write("mgRecoX");
  
  c1->cd(2);
  c1->cd(2)->SetGridx();
  c1->cd(2)->SetGridy();
  mgReco_XYZ[1]->GetXaxis()->SetLabelSize(0.05);
  mgReco_XYZ[1]->GetYaxis()->SetLabelSize(0.06);
  mgReco_XYZ[1]->GetXaxis()->SetTitleSize(0.05);
  mgReco_XYZ[1]->GetYaxis()->SetTitleSize(0.06);
  mgReco_XYZ[1]->GetYaxis()->SetTitleOffset(0.7);
  mgReco_XYZ[1]->Draw("AP");
  mgReco_XYZ[1]->SetTitle(";Unixtime (s);Reconstructed Y (m);");
  OutputFile->cd();
  mgReco_XYZ[1]->Write("mgRecoY");
  
  c1->cd(3);
  c1->cd(3)->SetGridx();
  c1->cd(3)->SetGridy();
  mgReco_XYZ[2]->GetXaxis()->SetLabelSize(0.05);
  mgReco_XYZ[2]->GetYaxis()->SetLabelSize(0.06);
  mgReco_XYZ[2]->GetXaxis()->SetTitleSize(0.05);
  mgReco_XYZ[2]->GetYaxis()->SetTitleSize(0.06);
  mgReco_XYZ[2]->GetYaxis()->SetTitleOffset(0.7);
  mgReco_XYZ[2]->Draw("AP");
  mgReco_XYZ[2]->SetTitle(";Unixtime (s);Reconstructed Z (m);");  
  OutputFile->cd();
  mgReco_XYZ[2]->Write("mgRecoZ");

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
  mgReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  mgReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.06);
  mgReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  mgReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.06);
  mgReco_ThPhR[0]->GetYaxis()->SetTitleOffset(0.7);
  mgReco_ThPhR[0]->Draw("AP");
  mgReco_ThPhR[0]->GetXaxis()->SetTitle("Unixtime (s)");
  mgReco_ThPhR[0]->GetYaxis()->SetTitle("Reconstructed Theta (^{o})");
  OutputFile->cd();
  mgReco_ThPhR[0]->Write("mgRecoTh");
  
  c2->cd(2);
  c2->cd(2)->SetGridx();
  c2->cd(2)->SetGridy();
  mgReco_ThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  mgReco_ThPhR[1]->GetYaxis()->SetLabelSize(0.06);
  mgReco_ThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  mgReco_ThPhR[1]->GetYaxis()->SetTitleSize(0.06);
  mgReco_ThPhR[1]->GetYaxis()->SetTitleOffset(0.7);
  mgReco_ThPhR[1]->Draw("AP");
  mgReco_ThPhR[1]->GetXaxis()->SetTitle("Unixtime (s)");
  mgReco_ThPhR[1]->GetYaxis()->SetTitle("Reconstructed Phi (^{o})");
  OutputFile->cd();
  mgReco_ThPhR[1]->Write("mgRecoPh");
  
  c2->cd(3);
  c2->cd(3)->SetGridx();
  c2->cd(3)->SetGridy();
  mgReco_ThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  mgReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.06);
  mgReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  mgReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.06);
  mgReco_ThPhR[2]->GetYaxis()->SetTitleOffset(0.7);
  mgReco_ThPhR[2]->Draw("AP");
  mgReco_ThPhR[2]->GetXaxis()->SetTitle("Unixtime (s)");
  mgReco_ThPhR[2]->GetYaxis()->SetTitle("Reconstructed Radial Displacement (^{o})");
  OutputFile->cd();
  mgReco_ThPhR[2]->Write("mgRecoR");

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
  hReco_dXYZ[0]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dXYZ[1]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dXYZ[2]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dXYZ[0]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dXYZ[1]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dXYZ[2]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dThPhR[0]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetTitleSize(0.06); 
  hReco_dThPhR[0]->GetYaxis()->SetTitleOffset(0.7);
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
  hReco_dThPhR[1]->GetYaxis()->SetTitleOffset(0.7);
  
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
  hReco_dThPhR[2]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR[2]->Draw();
  hReco_dThPhR[2]->SetTitle(";Reconstructed Displacement True - Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_dThPhR[2]->Write("hReco_dR");

  c4->cd(2);
  c4->cd(2)->SetGridx();
  c4->cd(2)->SetGridy();
  //c4->cd(2)->SetLogy();
  hReco_dThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[0]->GetYaxis()->SetTitleSize(0.06); 
  hReco_dThPhR[0]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR[0]->Draw();
  hReco_dThPhR[0]->SetTitle(";Reconstructed Theta True - Reco (^{o});No. of Events;");
  
  c4->cd(4);
  c4->cd(4)->SetGridx();
  c4->cd(4)->SetGridy();
  //c4->cd(4)->SetLogy();
  hReco_dThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[1]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[1]->GetYaxis()->SetTitleSize(0.06);
  hReco_dThPhR[1]->GetYaxis()->SetTitleOffset(0.7); 
  hReco_dThPhR[1]->Draw();
  hReco_dThPhR[1]->SetTitle(";Reconstructed Phi True - Reco (^{o});No. of Events;");
 
  c4->cd(6);
  c4->cd(6)->SetGridx();
  c4->cd(6)->SetGridy();
  //c4->cd(6)->SetLogy();
  hReco_dThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR[2]->GetYaxis()->SetLabelSize(0.06);
  hReco_dThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR[2]->GetYaxis()->SetTitleSize(0.06);
  hReco_dThPhR[2]->GetYaxis()->SetTitleOffset(0.7);
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

  // TCanvas *c4b = new TCanvas("c4b","c4b",1800,1800);
  // c4b->Divide(1,3);
  // c4b->cd(1);
  // c4b->cd(1)->SetGridx();
  // c4b->cd(1)->SetGridy();
  // c4b->cd(1)->SetLogy();
  // hReco_dThPhR_fR[0]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dThPhR_fR[0]->GetYaxis()->SetLabelSize(0.06);
  // hReco_dThPhR_fR[0]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dThPhR_fR[0]->GetYaxis()->SetTitleSize(0.06);
  // hReco_dThPhR_fR[0]->Draw();
  // hReco_dThPhR_fR[0]->SetTitle(";Reconstructed Theta True - Reco (Fixed R) (^{o});No. of Events;");
  // OutputFile->cd();
  // hReco_dThPhR_fR[0]->Write("hReco_dTh_fR");
  
  // c4b->cd(2);
  // c4b->cd(2)->SetGridx();
  // c4b->cd(2)->SetGridy();
  // c4b->cd(2)->SetLogy();
  // hReco_dThPhR_fR[1]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dThPhR_fR[1]->GetYaxis()->SetLabelSize(0.06);
  // hReco_dThPhR_fR[1]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dThPhR_fR[1]->GetYaxis()->SetTitleSize(0.06);
  // hReco_dThPhR_fR[1]->Draw();
  // hReco_dThPhR_fR[1]->SetTitle(";Reconstructed Phi True - Reco (Fixed R)(^{o});No. of Events;");
  // OutputFile->cd();
  // hReco_dThPhR_fR[1]->Write("hReco_dPh_fR");

  // c4b->cd(3);
  // c4b->cd(3)->SetGridx();
  // c4b->cd(3)->SetGridy();
  // c4b->cd(3)->SetLogy();
  // hReco_dThPhR_fR[2]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dThPhR_fR[2]->GetYaxis()->SetLabelSize(0.06); 
  // hReco_dThPhR_fR[2]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dThPhR_fR[2]->GetYaxis()->SetTitleSize(0.06);
  // hReco_dThPhR_fR[2]->Draw();
  // hReco_dThPhR_fR[2]->SetTitle(";Reconstructed Displacement True - Reco (Fixed R)(m);No. of Events;");
  // OutputFile->cd();
  // hReco_dThPhR_fR[2]->Write("hReco_dR_fR");

  // PlotFileName="Run";
  // PlotFileName+=Run;
  // PlotFileName+="RecodThPhR_fR";
  // PlotFileName+=".png";
  // c4b->SaveAs(PlotFileName);
  // stringL=PlotFileName.Length();
  // PlotFileName.Replace(stringL-4,4,".pdf");
  // c4b->SaveAs(PlotFileName);
  
  TCanvas *c5 = new TCanvas("c5","c5",1800,1800);
  c5->cd(1);
  c5->cd(1)->SetGridx();
  c5->cd(1)->SetGridy();
  mgDuration->GetXaxis()->SetLabelSize(0.05);
  mgDuration->GetYaxis()->SetLabelSize(0.06);
  mgDuration->GetXaxis()->SetTitleSize(0.05);
  mgDuration->GetYaxis()->SetTitleSize(0.06);
  mgDuration->GetYaxis()->SetTitleOffset(0.7);
  mgDuration->Draw("AP");
  mgDuration->GetXaxis()->SetTitle("Run Unixtime (s)");
  mgDuration->GetYaxis()->SetTitle("Time Taken (s)");
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
  grIterationsReco->GetYaxis()->SetLabelSize(0.06);
  grIterationsReco->GetXaxis()->SetTitleSize(0.05);
  grIterationsReco->GetYaxis()->SetTitleSize(0.06);
  grIterationsReco->GetYaxis()->SetTitleOffset(0.7);
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

}
