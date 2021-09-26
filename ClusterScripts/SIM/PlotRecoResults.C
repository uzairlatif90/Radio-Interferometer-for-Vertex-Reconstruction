void PlotRecoResults(){

  TString Run="Sim";

  double eventNum;
  double DurationTotal;
  double DurationReconstruction;
  double DurationInitialCondition;
  double FinalMinValue;
  int Iterations;
  double FinalTxCor_XYZ[3];
  double FinalTxCor_ThPhR[3];
  double InitialTxCor_XYZ[3];
  double InitialTxCor_ThPhR[3];
  double TrueTxCor_XYZ[3];
  double TrueTxCor_ThPhR[3];
 
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
  RecoTree->SetBranchAddress("DurationTotal",&DurationTotal);
  RecoTree->SetBranchAddress("DurationReconstruction",&DurationReconstruction);
  RecoTree->SetBranchAddress("DurationInitialCondition",&DurationInitialCondition);
  RecoTree->SetBranchAddress("FinalMinValue",&FinalMinValue);
  RecoTree->SetBranchAddress("Iterations",&Iterations);

  RecoTree->SetBranchAddress("FinalTxCor_XYZ",FinalTxCor_XYZ);
  RecoTree->SetBranchAddress("FinalTxCor_ThPhR",FinalTxCor_ThPhR);
  RecoTree->SetBranchAddress("InitialTxCor_XYZ",InitialTxCor_XYZ);
  RecoTree->SetBranchAddress("InitialTxCor_ThPhR",InitialTxCor_ThPhR);  
  RecoTree->SetBranchAddress("TrueTxCor_XYZ",TrueTxCor_XYZ);
  RecoTree->SetBranchAddress("TrueTxCor_ThPhR",TrueTxCor_ThPhR);  

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

  TH1D * hReco_dThPhR_Ini[3];
  hReco_dThPhR_Ini[0]=new TH1D("","",200,-200,200);
  hReco_dThPhR_Ini[1]=new TH1D("","",200,-200,200);  
  hReco_dThPhR_Ini[2]=new TH1D("","",200,-500,500);

  TH1D * hReco_dThPhR_True[3];
  hReco_dThPhR_True[0]=new TH1D("","",200,-200,200);
  hReco_dThPhR_True[1]=new TH1D("","",200,-200,200);  
  hReco_dThPhR_True[2]=new TH1D("","",200,-500,500);

  TH1D * hReco_ThPhR[3];
  hReco_ThPhR[0]=new TH1D("","",200,0,90);
  hReco_ThPhR[1]=new TH1D("","",200,-180,180);  
  hReco_ThPhR[2]=new TH1D("","",200,0,1001);

  TH2D * hDurationTotal[3];
  hDurationTotal[0]=new TH2D("","",100,90,180,100,0,300);
  hDurationTotal[1]=new TH2D("","",100,-180,180,100,0,300);
  hDurationTotal[2]=new TH2D("","",100,0,1010,100,0,300);

  TH2D * hDurationReco[3];
  hDurationReco[0]=new TH2D("","",100,90,180,100,0,300);
  hDurationReco[1]=new TH2D("","",100,-180,180,100,0,300);
  hDurationReco[2]=new TH2D("","",100,0,1010,100,0,300);
  
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
    if(FinalTxCor_XYZ[0]!=0 && FinalTxCor_XYZ[1]!=0 && TrueTxCor_ThPhR[0]<90){
      hIterations->Fill(Iterations);
      IterationsReco.push_back(Iterations);
      DurationReco.push_back(DurationReconstruction/1000);
      DurationTot.push_back(DurationTotal/1000);
      DurationIniCon.push_back(DurationInitialCondition/1000);
      UnixTimeReco.push_back(eventNum);

      hDurationTotal[0]->Fill(TrueTxCor_ThPhR[0],DurationTotal/1000);
      hDurationTotal[1]->Fill(TrueTxCor_ThPhR[1],DurationTotal/1000);
      hDurationTotal[2]->Fill(TrueTxCor_ThPhR[2],DurationTotal/1000);

      hDurationReco[0]->Fill(TrueTxCor_ThPhR[0],DurationReconstruction/1000);
      hDurationReco[1]->Fill(TrueTxCor_ThPhR[1],DurationReconstruction/1000);
      hDurationReco[2]->Fill(TrueTxCor_ThPhR[2],DurationReconstruction/1000);
      
      //cout<<setprecision(10)<<ievt<<" "<<unixTime<<" "<<firstUnixTime<<" "<<unixTime-firstUnixTime<<endl;
    
      Reco_XYZ[0].push_back(FinalTxCor_XYZ[0]);
      Reco_XYZ[1].push_back(FinalTxCor_XYZ[1]);
      Reco_XYZ[2].push_back(FinalTxCor_XYZ[2]);
      hReco_dXYZ[0]->Fill(InitialTxCor_XYZ[0]-FinalTxCor_XYZ[0]);
      hReco_dXYZ[1]->Fill(InitialTxCor_XYZ[1]-FinalTxCor_XYZ[1]);
      hReco_dXYZ[2]->Fill(InitialTxCor_XYZ[2]-FinalTxCor_XYZ[2]);

      Reco_ThPhR[0].push_back(FinalTxCor_ThPhR[0]);
      Reco_ThPhR[1].push_back(FinalTxCor_ThPhR[1]);
      Reco_ThPhR[2].push_back(FinalTxCor_ThPhR[2]);
      
      hReco_dThPhR_True[0]->Fill(TrueTxCor_ThPhR[0]-FinalTxCor_ThPhR[0]);
      hReco_dThPhR_True[1]->Fill(TrueTxCor_ThPhR[1]-FinalTxCor_ThPhR[1]);
      hReco_dThPhR_True[2]->Fill(TrueTxCor_ThPhR[2]-FinalTxCor_ThPhR[2]);

      hReco_dThPhR_Ini[0]->Fill(-InitialTxCor_ThPhR[0]+TrueTxCor_ThPhR[0]);
      hReco_dThPhR_Ini[1]->Fill(-InitialTxCor_ThPhR[1]+TrueTxCor_ThPhR[1]);
      hReco_dThPhR_Ini[2]->Fill(-InitialTxCor_ThPhR[2]+TrueTxCor_ThPhR[2]);
     
      
      True_XYZ[0].push_back(TrueTxCor_XYZ[0]);
      True_XYZ[1].push_back(TrueTxCor_XYZ[1]);
      True_XYZ[2].push_back(TrueTxCor_XYZ[2]);

      True_ThPhR[0].push_back(TrueTxCor_ThPhR[0]);
      True_ThPhR[1].push_back(TrueTxCor_ThPhR[1]);
      True_ThPhR[2].push_back(TrueTxCor_ThPhR[2]);

      if(fabs(TrueTxCor_ThPhR[0]-FinalTxCor_ThPhR[0])>80){
      	//cout<<TrueTxCor_ThPhR[0]<<" "<<TrueTxCor_ThPhR[1]<<" "<<TrueTxCor_ThPhR[2]<<endl;
	hReco_ThPhR[0]->Fill(TrueTxCor_ThPhR[0]);
	hReco_ThPhR[1]->Fill(TrueTxCor_ThPhR[1]);
	hReco_ThPhR[2]->Fill(TrueTxCor_ThPhR[2]);
      }
    }
  }////event loop

  // TMultiGraph * mgReco_XYZ[3];
  // mgReco_XYZ[0]=new TMultiGraph();
  // mgReco_XYZ[1]=new TMultiGraph();
  // mgReco_XYZ[2]=new TMultiGraph();
  // TGraph *grReco_XYZ[3];
  // TGraph *grTrue_XYZ[3];

  // TMultiGraph * mgReco_ThPhR[3];
  // mgReco_ThPhR[0]=new TMultiGraph();
  // mgReco_ThPhR[1]=new TMultiGraph();
  // mgReco_ThPhR[2]=new TMultiGraph();
  // TGraph *grReco_ThPhR[3];
  // TGraph *grTrue_ThPhR[3];
  
  // for(Int_t ixyz=0; ixyz<3; ixyz++){
  //   grReco_XYZ[ixyz]=new TGraph(Reco_XYZ[ixyz].size(),UnixTimeReco.data(),Reco_XYZ[ixyz].data());
  //   grReco_XYZ[ixyz]->SetMarkerStyle(20);
  //   grReco_XYZ[ixyz]->SetMarkerColor(2+ixyz);

  //   grReco_ThPhR[ixyz]=new TGraph(Reco_ThPhR[ixyz].size(),UnixTimeReco.data(),Reco_ThPhR[ixyz].data());
  //   grReco_ThPhR[ixyz]->SetMarkerStyle(20);
  //   grReco_ThPhR[ixyz]->SetMarkerColor(2+ixyz);

  //   grTrue_XYZ[ixyz]=new TGraph(True_XYZ[ixyz].size(),UnixTimeReco.data(),True_XYZ[ixyz].data());
  //   grTrue_XYZ[ixyz]->SetMarkerStyle(43);
  //   grTrue_XYZ[ixyz]->SetMarkerColor(1);

  //   grTrue_ThPhR[ixyz]=new TGraph(True_ThPhR[ixyz].size(),UnixTimeReco.data(),True_ThPhR[ixyz].data());
  //   grTrue_ThPhR[ixyz]->SetMarkerStyle(43);
  //   grTrue_ThPhR[ixyz]->SetMarkerColor(1);

  //   mgReco_XYZ[ixyz]->Add(grReco_XYZ[ixyz]);
  //   mgReco_XYZ[ixyz]->Add(grTrue_XYZ[ixyz]);

  //   mgReco_ThPhR[ixyz]->Add(grReco_ThPhR[ixyz]);
  //   mgReco_ThPhR[ixyz]->Add(grTrue_ThPhR[ixyz]);
    
  // }

  TGraph *grDurationReco[3];
  TGraph *grDurationTot[3];
  TGraph *grDurationIniCon[3];
  TMultiGraph *mgDuration[3];
  
  for(Int_t ixyz=0; ixyz<3; ixyz++){
    grDurationReco[ixyz]=new TGraph(DurationReco.size(),True_ThPhR[ixyz].data(),DurationReco.data());
    grDurationReco[ixyz]->SetMarkerStyle(20);
    grDurationReco[ixyz]->SetMarkerColor(kRed);
    grDurationReco[ixyz]->SetMarkerSize(1);
    grDurationReco[ixyz]->SetTitle("Reconstruction Time (s)");  

    grDurationTot[ixyz]=new TGraph(DurationTot.size(),True_ThPhR[ixyz].data(),DurationTot.data());
    grDurationTot[ixyz]->SetMarkerStyle(20);
    grDurationTot[ixyz]->SetMarkerColor(kBlue);
    grDurationTot[ixyz]->SetMarkerSize(3);
    grDurationTot[ixyz]->SetTitle("Total Time (s)");
  
    grDurationIniCon[ixyz]=new TGraph(DurationIniCon.size(),True_ThPhR[ixyz].data(),DurationIniCon.data());
    grDurationIniCon[ixyz]->SetMarkerStyle(20);
    grDurationIniCon[ixyz]->SetMarkerColor(kGreen);
    grDurationIniCon[ixyz]->SetMarkerSize(ixyz+1);
    grDurationIniCon[ixyz]->SetMarkerSize(2);
    grDurationIniCon[ixyz]->SetTitle("Initial Condition Time (s)");  

    mgDuration[ixyz]=new TMultiGraph();
    mgDuration[ixyz]->Add(grDurationTot[ixyz]);
    mgDuration[ixyz]->Add(grDurationIniCon[ixyz]);
    mgDuration[ixyz]->Add(grDurationReco[ixyz]);
  }

  TGraph *grIterationsReco=new TGraph(IterationsReco.size(),UnixTimeReco.data(),IterationsReco.data());
  grIterationsReco->SetMarkerStyle(20);
  grIterationsReco->SetMarkerColor(2);

  
  // TCanvas *c1 = new TCanvas("c1","c1",1800,1800);
  // c1->Divide(1,3);
  // c1->cd(1);
  // c1->cd(1)->SetGridx();
  // c1->cd(1)->SetGridy();
  // mgReco_XYZ[0]->GetXaxis()->SetLabelSize(0.05);
  // mgReco_XYZ[0]->GetYaxis()->SetLabelSize(0.05);
  // mgReco_XYZ[0]->GetXaxis()->SetTitleSize(0.05);
  // mgReco_XYZ[0]->GetYaxis()->SetTitleSize(0.05);
  // mgReco_XYZ[0]->Draw("AP");
  // mgReco_XYZ[0]->SetTitle(";Unixtime (s);Reconstructed X (m);");
  // OutputFile->cd();
  // mgReco_XYZ[0]->Write("mgRecoX");
  
  // c1->cd(2);
  // c1->cd(2)->SetGridx();
  // c1->cd(2)->SetGridy();
  // mgReco_XYZ[1]->GetXaxis()->SetLabelSize(0.05);
  // mgReco_XYZ[1]->GetYaxis()->SetLabelSize(0.05);
  // mgReco_XYZ[1]->GetXaxis()->SetTitleSize(0.05);
  // mgReco_XYZ[1]->GetYaxis()->SetTitleSize(0.05);
  // mgReco_XYZ[1]->Draw("AP");
  // mgReco_XYZ[1]->SetTitle(";Unixtime (s);Reconstructed Y (m);");
  // OutputFile->cd();
  // mgReco_XYZ[1]->Write("mgRecoY");
  
  // c1->cd(3);
  // c1->cd(3)->SetGridx();
  // c1->cd(3)->SetGridy();
  // mgReco_XYZ[2]->GetXaxis()->SetLabelSize(0.05);
  // mgReco_XYZ[2]->GetYaxis()->SetLabelSize(0.05);
  // mgReco_XYZ[2]->GetXaxis()->SetTitleSize(0.05);
  // mgReco_XYZ[2]->GetYaxis()->SetTitleSize(0.05);
  // mgReco_XYZ[2]->Draw("AP");
  // mgReco_XYZ[2]->SetTitle(";Unixtime (s);Reconstructed Z (m);");  
  // OutputFile->cd();
  // mgReco_XYZ[2]->Write("mgRecoZ");

  TString PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecoXYZ";
  PlotFileName+=".png";
  //c1->SaveAs(PlotFileName);
  int stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  //c1->SaveAs(PlotFileName);
  
  // TCanvas *c2 = new TCanvas("c2","c2",1800,1800);
  // c2->Divide(1,3);
  // c2->cd(1);
  // c2->cd(1)->SetGridx();
  // c2->cd(1)->SetGridy();
  // mgReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.1);
  // mgReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.1);
  // mgReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.1);
  // mgReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.1);

  // mgReco_ThPhR[0]->Draw("AP");
  // mgReco_ThPhR[0]->SetTitle(";Unixtime (s);Reconstructed Theta (^{o});");
  // OutputFile->cd();
  // mgReco_ThPhR[0]->Write("mgRecoTh");
  
  // c2->cd(2);
  // c2->cd(2)->SetGridx();
  // c2->cd(2)->SetGridy();
  // mgReco_ThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  // mgReco_ThPhR[1]->GetYaxis()->SetLabelSize(0.05);
  // mgReco_ThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  // mgReco_ThPhR[1]->GetYaxis()->SetTitleSize(0.05);
  // mgReco_ThPhR[1]->Draw("AP");
  // mgReco_ThPhR[1]->SetTitle(";Unixtime (s);Reconstructed Phi (^{o});");
  // OutputFile->cd();
  // mgReco_ThPhR[1]->Write("mgRecoPh");
  
  // c2->cd(3);
  // c2->cd(3)->SetGridx();
  // c2->cd(3)->SetGridy();
  // mgReco_ThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  // mgReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.05);
  // mgReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  // mgReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.05);
  // mgReco_ThPhR[2]->Draw("AP");
  // mgReco_ThPhR[2]->SetTitle(";Unixtime (s);Reconstructed Displacement (m);");  
  // OutputFile->cd();
  // mgReco_ThPhR[2]->Write("mgRecoR");

  // PlotFileName="Run";
  // PlotFileName+=Run;
  // PlotFileName+="RecoThPhR";
  // PlotFileName+=".png";
  // c2->SaveAs(PlotFileName);
  // stringL=PlotFileName.Length();
  // PlotFileName.Replace(stringL-4,4,".pdf");
  // // c2->SaveAs(PlotFileName);
  
  // TCanvas *c3 = new TCanvas("c3","c3",1800,1800);
  // c3->Divide(1,3);
  // c3->cd(1);
  // c3->cd(1)->SetGridx();
  // c3->cd(1)->SetGridy();
  // c3->cd(1)->SetLogy();
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
  // c3->cd(2)->SetLogy();
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
  // c3->cd(3)->SetLogy();
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

  TCanvas *c4 = new TCanvas("c4","c4",1800,1800);
  c4->Divide(2,3);
  c4->cd(1);
  c4->cd(1)->SetGridx();
  c4->cd(1)->SetGridy();
  c4->cd(1)->SetLogy();
  hReco_dThPhR_Ini[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_Ini[0]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_Ini[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_Ini[0]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_Ini[0]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_Ini[0]->Draw();
  hReco_dThPhR_Ini[0]->SetTitle(";Theta True - Initial (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_Ini[0]->Write("hReco_dTh_Initial");
  
  c4->cd(3);
  c4->cd(3)->SetGridx();
  c4->cd(3)->SetGridy();
  c4->cd(3)->SetLogy();
  hReco_dThPhR_Ini[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_Ini[1]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_Ini[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_Ini[1]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_Ini[1]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_Ini[1]->Draw();
  hReco_dThPhR_Ini[1]->SetTitle(";Phi True - Initial (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_Ini[1]->Write("hReco_dPh_Initial");

  c4->cd(5);
  c4->cd(5)->SetGridx();
  c4->cd(5)->SetGridy();
  c4->cd(5)->SetLogy();
  hReco_dThPhR_Ini[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_Ini[2]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_Ini[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_Ini[2]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_Ini[2]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_Ini[2]->Draw();
  hReco_dThPhR_Ini[2]->SetTitle(";Displacement True - Initial (m);No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_Ini[2]->Write("hReco_dR_Initial");

  c4->cd(2);
  c4->cd(2)->SetGridx();
  c4->cd(2)->SetGridy();
  //c4->cd(2)->SetLogy();
  hReco_dThPhR_Ini[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_Ini[0]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_Ini[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_Ini[0]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_Ini[0]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_Ini[0]->Draw();
  hReco_dThPhR_Ini[0]->SetTitle(";Theta True - Initial (^{o});No. of Events;");
  
  c4->cd(4);
  c4->cd(4)->SetGridx();
  c4->cd(4)->SetGridy();
  //c4->cd(4)->SetLogy();
  hReco_dThPhR_Ini[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_Ini[1]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_Ini[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_Ini[1]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_Ini[1]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_Ini[1]->Draw();
  hReco_dThPhR_Ini[1]->SetTitle(";Phi True - Initial (^{o});No. of Events;");
 
  c4->cd(6);
  c4->cd(6)->SetGridx();
  c4->cd(6)->SetGridy();
  //c4->cd(6)->SetLogy();
  hReco_dThPhR_Ini[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_Ini[2]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_Ini[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_Ini[2]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_Ini[2]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_Ini[2]->Draw();
  hReco_dThPhR_Ini[2]->SetTitle(";Displacement True - Initial (m);No. of Events;");
 
  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecodThPhR_Initial";
  PlotFileName+=".png";
  c4->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c4->SaveAs(PlotFileName);

  TCanvas *c4b = new TCanvas("c4b","c4b",1800,1800);
  c4b->Divide(2,3);
  c4b->cd(1);
  c4b->cd(1)->SetGridx();
  c4b->cd(1)->SetGridy();
  c4b->cd(1)->SetLogy();
  hReco_dThPhR_True[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_True[0]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_True[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_True[0]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_True[0]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_True[0]->Draw();
  hReco_dThPhR_True[0]->SetTitle(";Reconstructed Theta True - Reco (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_True[0]->Write("hReco_dTh_True");
  
  c4b->cd(3);
  c4b->cd(3)->SetGridx();
  c4b->cd(3)->SetGridy();
  c4b->cd(3)->SetLogy();
  hReco_dThPhR_True[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_True[1]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_True[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_True[1]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_True[1]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_True[1]->Draw();
  hReco_dThPhR_True[1]->SetTitle(";Reconstructed Phi True - Reco (^{o});No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_True[1]->Write("hReco_dPh_True");

  c4b->cd(5);
  c4b->cd(5)->SetGridx();
  c4b->cd(5)->SetGridy();
  c4b->cd(5)->SetLogy();
  hReco_dThPhR_True[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_True[2]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_True[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_True[2]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_True[2]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_True[2]->Draw();
  hReco_dThPhR_True[2]->SetTitle(";Reconstructed Displacement True - Reco (m);No. of Events;");
  OutputFile->cd();
  hReco_dThPhR_True[2]->Write("hReco_dR_True");

  c4b->cd(2);
  c4b->cd(2)->SetGridx();
  c4b->cd(2)->SetGridy();
  //c4b->cd(2)->SetLogy();
  hReco_dThPhR_True[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_True[0]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_True[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_True[0]->GetYaxis()->SetTitleSize(0.07); 
 hReco_dThPhR_True[0]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_True[0]->Draw();
  hReco_dThPhR_True[0]->SetTitle(";Reconstructed Theta True - Reco (^{o});No. of Events;");
  
  c4b->cd(4);
  c4b->cd(4)->SetGridx();
  c4b->cd(4)->SetGridy();
  //c4b->cd(4)->SetLogy();
  hReco_dThPhR_True[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_True[1]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_True[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_True[1]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_True[1]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_True[1]->Draw();
  hReco_dThPhR_True[1]->SetTitle(";Reconstructed Phi True - Reco (^{o});No. of Events;");
 
  c4b->cd(6);
  c4b->cd(6)->SetGridx();
  c4b->cd(6)->SetGridy();
  //c4b->cd(6)->SetLogy();
  hReco_dThPhR_True[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_dThPhR_True[2]->GetYaxis()->SetLabelSize(0.07);
  hReco_dThPhR_True[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_dThPhR_True[2]->GetYaxis()->SetTitleSize(0.07);
  hReco_dThPhR_True[2]->GetYaxis()->SetTitleOffset(0.7);
  hReco_dThPhR_True[2]->Draw();
  hReco_dThPhR_True[2]->SetTitle(";Reconstructed Displacement True - Reco (m);No. of Events;");
 
  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecodThPhR_True";
  PlotFileName+=".png";
  c4b->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c4b->SaveAs(PlotFileName);
  
  TCanvas *c5 = new TCanvas("c5","c5",1800,1800);
  c5->Divide(1,3);
  for(int ixyz=0;ixyz<3;ixyz++){
    c5->cd(1+ixyz);
    c5->cd(1+ixyz)->SetGridx();
    c5->cd(1+ixyz)->SetGridy();
    mgDuration[ixyz]->GetXaxis()->SetLabelSize(0.05);
    mgDuration[ixyz]->GetYaxis()->SetLabelSize(0.07);
    mgDuration[ixyz]->GetXaxis()->SetTitleSize(0.05);
    mgDuration[ixyz]->GetYaxis()->SetTitleSize(0.07);
    mgDuration[ixyz]->GetYaxis()->SetTitleOffset(0.7);
    mgDuration[ixyz]->GetXaxis()->SetNdivisions(20);

    // mgDuration[ixyz]->GetHistogram()->SetLabelSize(0.09);
    // mgDuration[ixyz]->GetHistogram()->SetLabelSize(0.09);
    // mgDuration[ixyz]->GetHistogram()->SetTitleSize(0.09);
    // mgDuration[ixyz]->GetHistogram()->SetTitleSize(0.09);
  
    mgDuration[ixyz]->Draw("AP");
    TString title="";
    if(ixyz==0){
      title="Theta (deg)";
    }
    if(ixyz==1){
      title="Phi (deg)";
    }
    if(ixyz==2){
      title="Radial Displacement (m)";
    }
    mgDuration[ixyz]->GetHistogram()->SetTitle(title);
    //mgDuration[ixyz]->GetXaxis()->SetTitle(title);
    mgDuration[ixyz]->GetYaxis()->SetTitle("Time Taken (s)");
    c5->cd(1+ixyz)->BuildLegend();
    OutputFile->cd();
    title+="mgDuration";
    mgDuration[ixyz]->Write(title);

    PlotFileName="Run";
    PlotFileName+=Run;
    PlotFileName+="Durations";
    PlotFileName+=".png";
    c5->SaveAs(PlotFileName);
    stringL=PlotFileName.Length();
    PlotFileName.Replace(stringL-4,4,".pdf");
    c5->SaveAs(PlotFileName);
  }
  
  TCanvas *c6 = new TCanvas("c6","c6",1800,1800);
  c6->cd(1);
  c6->cd(1)->SetGridx();
  c6->cd(1)->SetGridy();
  // grIterationsReco->GetXaxis()->SetLabelSize(0.05);
  // grIterationsReco->GetYaxis()->SetLabelSize(0.06);
  // grIterationsReco->GetXaxis()->SetTitleSize(0.05);
  // grIterationsReco->GetYaxis()->SetTitleSize(0.06);
  // grIterationsReco->Draw("AP");
  // grIterationsReco->SetTitle(";Run Unixtime (s);Iterations taken by Minimizer to do Reco.;");  
  // OutputFile->cd();
  // grIterationsReco->Write("grRecoIterations");

  hIterations->GetXaxis()->SetLabelSize(0.05);
  hIterations->GetYaxis()->SetLabelSize(0.05);
  hIterations->GetXaxis()->SetTitleSize(0.05);
  hIterations->GetYaxis()->SetTitleSize(0.05);
  hIterations->GetYaxis()->SetTitleOffset(0.7);
  hIterations->Draw("");
  hIterations->SetTitle(";Iterations taken by Minimizer to do Reco.; No. of Events");  
  OutputFile->cd();
  hIterations->Write("hIterations");

  
  PlotFileName="Run";
  PlotFileName+=Run;
  PlotFileName+="RecoIterations";
  PlotFileName+=".png";
  c6->SaveAs(PlotFileName);
  stringL=PlotFileName.Length();
  PlotFileName.Replace(stringL-4,4,".pdf");
  c6->SaveAs(PlotFileName);

  TCanvas *c7 = new TCanvas("c7","c7",1800,1800);
  c7->Divide(1,3);
  for(int ixyz=0;ixyz<3;ixyz++){
    c7->cd(1+ixyz);
    c7->cd(1+ixyz)->SetGridx();
    c7->cd(1+ixyz)->SetGridy();
    hDurationTotal[ixyz]->GetXaxis()->SetLabelSize(0.05);
    hDurationTotal[ixyz]->GetYaxis()->SetLabelSize(0.07);
    hDurationTotal[ixyz]->GetXaxis()->SetTitleSize(0.05);
    hDurationTotal[ixyz]->GetYaxis()->SetTitleSize(0.07);
    hDurationTotal[ixyz]->GetYaxis()->SetTitleOffset(0.7);
   
    hDurationTotal[ixyz]->Draw("colz");
    TString title="";
    if(ixyz==0){
      title="Theta (deg)";
    }
    if(ixyz==1){
      title="Phi (deg)";
    }
    if(ixyz==2){
      title="Radial Displacement (m)";
    }
    //hDurationTotal[ixyz]->GetHistogram()->SetTitle(title);
    hDurationTotal[ixyz]->GetXaxis()->SetTitle(title);
    hDurationTotal[ixyz]->GetYaxis()->SetTitle("Total Time Taken (s)");
  
    OutputFile->cd();
    title+="hDurationTotal";
    hDurationTotal[ixyz]->Write(title);

    PlotFileName="Run";
    PlotFileName+=Run;
    PlotFileName+="hDurationTotal";
    PlotFileName+=".png";
    c7->SaveAs(PlotFileName);
    stringL=PlotFileName.Length();
    PlotFileName.Replace(stringL-4,4,".pdf");
    c7->SaveAs(PlotFileName);
  }

  TCanvas *c8 = new TCanvas("c8","c8",1800,1800);
  c8->Divide(1,3);
  for(int ixyz=0;ixyz<3;ixyz++){
    c8->cd(1+ixyz);
    c8->cd(1+ixyz)->SetGridx();
    c8->cd(1+ixyz)->SetGridy();
    hDurationReco[ixyz]->GetXaxis()->SetLabelSize(0.05);
    hDurationReco[ixyz]->GetYaxis()->SetLabelSize(0.07);
    hDurationReco[ixyz]->GetXaxis()->SetTitleSize(0.05);
    hDurationReco[ixyz]->GetYaxis()->SetTitleSize(0.07);
    hDurationReco[ixyz]->GetYaxis()->SetTitleOffset(0.7);
   
    hDurationReco[ixyz]->Draw("colz");
    TString title="";
    if(ixyz==0){
      title="Theta (deg)";
    }
    if(ixyz==1){
      title="Phi (deg)";
    }
    if(ixyz==2){
      title="Radial Displacement (m)";
    }
    //hDurationReco[ixyz]->GetHistogram()->SetTitle(title);
    hDurationReco[ixyz]->GetXaxis()->SetTitle(title);
    hDurationReco[ixyz]->GetYaxis()->SetTitle("Reco Time Taken (s)");
  
    OutputFile->cd();
    title+="hDurationReco";
    hDurationReco[ixyz]->Write(title);

    PlotFileName="Run";
    PlotFileName+=Run;
    PlotFileName+="hDurationReco";
    PlotFileName+=".png";
    c8->SaveAs(PlotFileName);
    stringL=PlotFileName.Length();
    PlotFileName.Replace(stringL-4,4,".pdf");
    c8->SaveAs(PlotFileName);
  }

  TCanvas *c9 = new TCanvas("c9","c9",1800,1800);
  c9->Divide(2,3);
  c9->cd(1);
  c9->cd(1)->SetGridx();
  c9->cd(1)->SetGridy();
  c9->cd(1)->SetLogy();
  hReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.07);
  hReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.07);
  hReco_ThPhR[0]->GetYaxis()->SetTitleOffset(0.7);
  hReco_ThPhR[0]->Draw();
  hReco_ThPhR[0]->SetTitle(";Theta True - Initial (^{o});No. of Events;");
  OutputFile->cd();
  hReco_ThPhR[0]->Write("hReco_Th");
  
  c9->cd(3);
  c9->cd(3)->SetGridx();
  c9->cd(3)->SetGridy();
  c9->cd(3)->SetLogy();
  hReco_ThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR[1]->GetYaxis()->SetLabelSize(0.07);
  hReco_ThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR[1]->GetYaxis()->SetTitleSize(0.07);
  hReco_ThPhR[1]->GetYaxis()->SetTitleOffset(0.7);
  hReco_ThPhR[1]->Draw();
  hReco_ThPhR[1]->SetTitle(";Phi True - Initial (^{o});No. of Events;");
  OutputFile->cd();
  hReco_ThPhR[1]->Write("hReco_dPh");

  c9->cd(5);
  c9->cd(5)->SetGridx();
  c9->cd(5)->SetGridy();
  c9->cd(5)->SetLogy();
  hReco_ThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  hReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.07);
  hReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  hReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.07);
  hReco_ThPhR[2]->GetYaxis()->SetTitleOffset(0.7);
  hReco_ThPhR[2]->Draw();
  hReco_ThPhR[2]->SetTitle(";Displacement True - Initial (m);No. of Events;");
  OutputFile->cd();
  hReco_ThPhR[2]->Write("hReco_dR");

  // c9->cd(2);
  // c9->cd(2)->SetGridx();
  // c9->cd(2)->SetGridy();
  // //c9->cd(2)->SetLogy();
  // hReco_ThPhR[0]->GetXaxis()->SetLabelSize(0.05);
  // hReco_ThPhR[0]->GetYaxis()->SetLabelSize(0.07);
  // hReco_ThPhR[0]->GetXaxis()->SetTitleSize(0.05);
  // hReco_ThPhR[0]->GetYaxis()->SetTitleSize(0.07);
  // hReco_ThPhR[0]->GetYaxis()->SetTitleOffset(0.7);
  // hReco_ThPhR[0]->Draw();
  // hReco_ThPhR[0]->SetTitle(";Theta True - Initial (^{o});No. of Events;");
  
  // c9->cd(4);
  // c9->cd(4)->SetGridx();
  // c9->cd(4)->SetGridy();
  // //c9->cd(4)->SetLogy();
  // hReco_ThPhR[1]->GetXaxis()->SetLabelSize(0.05);
  // hReco_ThPhR[1]->GetYaxis()->SetLabelSize(0.07);
  // hReco_ThPhR[1]->GetXaxis()->SetTitleSize(0.05);
  // hReco_ThPhR[1]->GetYaxis()->SetTitleSize(0.07);
  // hReco_ThPhR[1]->GetYaxis()->SetTitleOffset(0.7);
  // hReco_ThPhR[1]->Draw();
  // hReco_ThPhR[1]->SetTitle(";Phi True - Initial (^{o});No. of Events;");
 
  // c9->cd(6);
  // c9->cd(6)->SetGridx();
  // c9->cd(6)->SetGridy();
  // //c9->cd(6)->SetLogy();
  // hReco_ThPhR[2]->GetXaxis()->SetLabelSize(0.05);
  // hReco_ThPhR[2]->GetYaxis()->SetLabelSize(0.07);
  // hReco_ThPhR[2]->GetXaxis()->SetTitleSize(0.05);
  // hReco_ThPhR[2]->GetYaxis()->SetTitleSize(0.07);
  // hReco_ThPhR[2]->GetYaxis()->SetTitleOffset(0.7);
  // hReco_ThPhR[2]->Draw();
  // hReco_ThPhR[2]->SetTitle(";Displacement True - Initial (m);No. of Events;");
 
  // PlotFileName="Run";
  // PlotFileName+=Run;
  // PlotFileName+="RecoThPhR";
  // PlotFileName+=".png";
  // c9->SaveAs(PlotFileName);
  // stringL=PlotFileName.Length();
  // PlotFileName.Replace(stringL-4,4,".pdf");
  // c9->SaveAs(PlotFileName);
  
}
