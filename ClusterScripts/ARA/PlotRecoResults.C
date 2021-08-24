void PlotRecoResults(){

  Int_t Run=9129;

  int eventNum;
  double unixTime;
  double firstUnixTime;
  double timeStamp;
  int runNum;
  double Duration;
  double FinalMinValue;
  int Iterations;
  bool isCalpulserTrig;
  bool isSoftwareTrig;
  double FinalTxCor_XYZ[3];
  double FinalTxCor_ThPhR[3];
  double InitialTxCor_XYZ[3];
  double InitialTxCor_ThPhR[3];
  double VoltageSNR[16]; 
  double VoltageSNRV[3]; 
  int VoltageSNRV_Ch[3];
  double VoltageSNRH[3]; 
  int VoltageSNRH_Ch[3];
  double CorScore[16];
  double PowerSNR[16];  
 
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
  RecoTree->SetBranchAddress("Duration",&Duration);
  RecoTree->SetBranchAddress("FinalMinValue",&FinalMinValue);
  RecoTree->SetBranchAddress("Iterations",&Iterations);
  RecoTree->SetBranchAddress("isCalpulserTrig",&isCalpulserTrig);
  RecoTree->SetBranchAddress("isSoftwareTrig",&isSoftwareTrig);
  RecoTree->SetBranchAddress("FinalTxCor_XYZ",FinalTxCor_XYZ);
  RecoTree->SetBranchAddress("FinalTxCor_ThPhR",FinalTxCor_ThPhR);
  RecoTree->SetBranchAddress("InitialTxCor_XYZ",InitialTxCor_XYZ);
  RecoTree->SetBranchAddress("InitialTxCor_ThPhR",InitialTxCor_ThPhR);  
  RecoTree->SetBranchAddress("VoltageSNR",VoltageSNR); 
  // RecoTree->SetBranchAddress("VoltageSNRV",VoltageSNRV);
  // RecoTree->SetBranchAddress("VoltageSNRV_Ch",VoltageSNRV_Ch);
  // RecoTree->SetBranchAddress("VoltageSNRH",VoltageSNRH);
  // RecoTree->SetBranchAddress("VoltageSNRH_Ch",VoltageSNRH_Ch);
  RecoTree->SetBranchAddress("CorScore",CorScore);
  RecoTree->SetBranchAddress("PowerSNR",PowerSNR);

  double numEntries=RecoTree->GetEntries();

  RecoTree->GetEntry(0);

  TH1D *hIterations=new TH1D("","",200,0,50);
  
  TH1D * hReco_dXYZ[3];
  hReco_dXYZ[0]=new TH1D("","",200,-500,+500);
  hReco_dXYZ[1]=new TH1D("","",200,-500,+500);
  hReco_dXYZ[2]=new TH1D("","",200,-500,+500);

  TH1D * hReco_dThPhR[3];
  hReco_dThPhR[0]=new TH1D("","",200,-5,5);
  hReco_dThPhR[1]=new TH1D("","",200,-5,5);  
  hReco_dThPhR[2]=new TH1D("","",200,-100,100);

  TH1D * hReco_ThPhR[3];
  hReco_ThPhR[0]=new TH1D("","",200,90,180);
  hReco_ThPhR[1]=new TH1D("","",200,-180,180);  
  hReco_ThPhR[2]=new TH1D("","",200,0,2000);
  
  const Int_t NumCutCh=4;
  Double_t CutCh[NumCutCh]={8,15,11,13};

  vector <double> IterationsReco;
  vector <double> DurationReco;
  vector <double> UnixTimeReco;
  vector <double> Reco_XYZ[3];
  vector <double> Reco_ThPhR[3];
  
  ///Start looping over the events
  for(Int_t ievt=0; ievt<numEntries; ievt++){
 
    RecoTree->GetEntry(ievt);

    if(InitialTxCor_XYZ[0]!=0 && InitialTxCor_XYZ[1]!=0 && InitialTxCor_XYZ[2]!=0 && FinalTxCor_ThPhR[2]<300){
    
    hIterations->Fill(Iterations);
    IterationsReco.push_back(Iterations);
    DurationReco.push_back(Duration);
    UnixTimeReco.push_back(unixTime-firstUnixTime);

    //cout<<setprecision(10)<<ievt<<" "<<unixTime<<" "<<firstUnixTime<<" "<<unixTime-firstUnixTime<<endl;
    
    Reco_XYZ[0].push_back(FinalTxCor_XYZ[0]);
    Reco_XYZ[1].push_back(FinalTxCor_XYZ[1]);
    Reco_XYZ[2].push_back(FinalTxCor_XYZ[2]);
    // hReco_dXYZ[0]->Fill(dX);
    // hReco_dXYZ[1]->Fill(dY);
    // hReco_dXYZ[2]->Fill(dZ);

    Reco_ThPhR[0].push_back(FinalTxCor_ThPhR[0]);
    Reco_ThPhR[1].push_back(FinalTxCor_ThPhR[1]);
    Reco_ThPhR[2].push_back(FinalTxCor_ThPhR[2]);

    hReco_ThPhR[0]->Fill(FinalTxCor_ThPhR[0]);
    hReco_ThPhR[1]->Fill(FinalTxCor_ThPhR[1]);
    hReco_ThPhR[2]->Fill(FinalTxCor_ThPhR[2]);

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
  grDurationReco->SetMarkerColor(2);

  TGraph *grIterationsReco=new TGraph(IterationsReco.size(),UnixTimeReco.data(),IterationsReco.data());
  grIterationsReco->SetMarkerStyle(20);
  grIterationsReco->SetMarkerColor(2);
  
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
  grReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  grReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.05);
  grReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  grReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.05);
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
  grReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.05);
  grReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  grReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.05);
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
  
  // TCanvas *c3 = new TCanvas("c3","c3",1800,1800);
  // c3->Divide(1,3);
  // c3->cd(1);
  // c3->cd(1)->SetGridx();
  // c3->cd(1)->SetGridy();
  // hReco_dXYZ[0]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dXYZ[0]->GetYaxis()->SetLabelSize(0.05);
  // hReco_dXYZ[0]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dXYZ[0]->GetYaxis()->SetTitleSize(0.05);
  // hReco_dXYZ[0]->Draw();
  // hReco_dXYZ[0]->SetTitle(";Reconstructed X True - Reco (m);No. of Events;");
  // OutputFile->cd();
  // hReco_dXYZ[0]->Write("hReco_dX");
  
  // c3->cd(2);
  // c3->cd(2)->SetGridx();
  // c3->cd(2)->SetGridy();
  // hReco_dXYZ[1]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dXYZ[1]->GetYaxis()->SetLabelSize(0.05);
  // hReco_dXYZ[1]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dXYZ[1]->GetYaxis()->SetTitleSize(0.05);
  // hReco_dXYZ[1]->Draw();
  // hReco_dXYZ[1]->SetTitle(";Reconstructed Y True - Reco (m);No. of Events;");
  // OutputFile->cd();
  // hReco_dXYZ[1]->Write("hReco_dY");

  // c3->cd(3);
  // c3->cd(3)->SetGridx();
  // c3->cd(3)->SetGridy();
  // hReco_dXYZ[2]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dXYZ[2]->GetYaxis()->SetLabelSize(0.05);
  // hReco_dXYZ[2]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dXYZ[2]->GetYaxis()->SetTitleSize(0.05);
  // hReco_dXYZ[2]->Draw();
  // hReco_dXYZ[2]->SetTitle(";Reconstructed Z True - Reco (m);No. of Events;");
  // OutputFile->cd();
  // hReco_dXYZ[2]->Write("hReco_dZ");

  // PlotFileName="Run";
  // PlotFileName+=Run;
  // PlotFileName+="RecodXYZ";
  // PlotFileName+=".png";
  // c3->SaveAs(PlotFileName);
  // stringL=PlotFileName.Length();
  // PlotFileName.Replace(stringL-4,4,".pdf");
  // c3->SaveAs(PlotFileName);

  // TCanvas *c4 = new TCanvas("c4","c4",1800,1800);
  // c4->Divide(1,3);
  // c4->cd(1);
  // c4->cd(1)->SetGridx();
  // c4->cd(1)->SetGridy();
  // hReco_dThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dThPhR[0]->GetYaxis()->SetLabelSize(0.05);
  // hReco_dThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dThPhR[0]->GetYaxis()->SetTitleSize(0.05);
  // hReco_dThPhR[0]->Draw();
  // hReco_dThPhR[0]->SetTitle(";Reconstructed Theta True - Reco (^{o});No. of Events;");
  // OutputFile->cd();
  // hReco_dThPhR[0]->Write("hReco_dTh");
  
  // c4->cd(2);
  // c4->cd(2)->SetGridx();
  // c4->cd(2)->SetGridy();
  // hReco_dThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dThPhR[1]->GetYaxis()->SetLabelSize(0.05);
  // hReco_dThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dThPhR[1]->GetYaxis()->SetTitleSize(0.05);
  // hReco_dThPhR[1]->Draw();
  // hReco_dThPhR[1]->SetTitle(";Reconstructed Phi True - Reco (^{o});No. of Events;");
  // OutputFile->cd();
  // hReco_dThPhR[1]->Write("hReco_dPh");

  // c4->cd(3);
  // c4->cd(3)->SetGridx();
  // c4->cd(3)->SetGridy();
  // hReco_dThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  // hReco_dThPhR[2]->GetYaxis()->SetLabelSize(0.05);
  // hReco_dThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  // hReco_dThPhR[2]->GetYaxis()->SetTitleSize(0.05);
  // hReco_dThPhR[2]->Draw();
  // hReco_dThPhR[2]->SetTitle(";Reconstructed Displacement True - Reco (m);No. of Events;");
  // OutputFile->cd();
  // hReco_dThPhR[2]->Write("hReco_dR");

  // PlotFileName="Run";
  // PlotFileName+=Run;
  // PlotFileName+="RecodThPhR";
  // PlotFileName+=".png";
  // c4->SaveAs(PlotFileName);
  // stringL=PlotFileName.Length();
  // PlotFileName.Replace(stringL-4,4,".pdf");
  // c4->SaveAs(PlotFileName);
  
  TCanvas *c5 = new TCanvas("c5","c5",1800,1800);
  c5->cd(1);
  c5->cd(1)->SetGridx();
  c5->cd(1)->SetGridy();
  grDurationReco->GetXaxis()->SetLabelSize(0.05);
  grDurationReco->GetYaxis()->SetLabelSize(0.05);
  grDurationReco->GetXaxis()->SetTitleSize(0.05);
  grDurationReco->GetYaxis()->SetTitleSize(0.05);
  grDurationReco->Draw("AP");
  grDurationReco->SetTitle(";Run Unixtime (s);Time Taken To Do Reconstruction (ms);");  
  OutputFile->cd();
  grDurationReco->Write("grRecoDuration");

  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecoDuration";
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
  hIterations->Draw();

  TCanvas *c8 = new TCanvas("c8","c8",1800,1800);
  c8->Divide(1,3);
  c8->cd(1);
  c8->cd(1)->SetGridx();
  c8->cd(1)->SetGridy();
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
  
}
