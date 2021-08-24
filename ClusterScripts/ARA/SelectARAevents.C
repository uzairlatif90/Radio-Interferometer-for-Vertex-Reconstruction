const Int_t MCH=16;

#include "FFTtools.h"
#include "/data/user/ulatif/Interferometer/Interferometer.cc"

Double_t TrueX=-456.721;
Double_t TrueY=-2353;

TGraph *gCPtemp[2][16];

void ReadCPTemp(){  

  {
    TString filename="/data/user/ulatif/SPICE_Inter/";
    filename+="CP_D6VPol_A2.root";
    TFile *f = TFile::Open(filename, "READ");
    for(int ich=0;ich<MCH;ich++){
      TString GetGraph="gr";
      GetGraph+=ich;
      f->GetObject(GetGraph, gCPtemp[0][ich]);
    }
    delete f;
  }
  {
    TString filename="/data/user/ulatif/SPICE_Inter/";
    filename+="CP_D6HPol_A2.root";
    TFile *f = TFile::Open(filename, "READ");
    for(int ich=0;ich<MCH;ich++){
      TString GetGraph="gr";
      GetGraph+=ich;
      f->GetObject(GetGraph, gCPtemp[1][ich]);
    }
    delete f;
  }
  
}

void SelectARAevents(Int_t StationId, TString newfileString){

  ReadCPTemp();
  //Int_t StationId=2;
  DeclareAntennaConfigARA(StationId);

  ////initialise the event pointer to take data from the event tree
  RawAtriStationEvent *rawAtriEvPtr=0;

  ////initialise the AraEventCalibrator class  
  AraEventCalibrator *calib = AraEventCalibrator::Instance();
  
  ///Open the ARA root file that contains the data
  //TFile *newfile=new TFile("/data/exp/SPICEcore/spiceCoreData2018-2019/L1/ARA02/1226/run_012576/event012576.root");
  TFile *newfile=new TFile(newfileString);
   
  ///Open the Tree in the file that contains data from all of the events
  TTree *ceventTree = (TTree*)newfile->Get("eventTree");

  Int_t runNumber;
  ///Set the Branch addresses
  ceventTree->SetBranchAddress("event",&rawAtriEvPtr);
  ceventTree->SetBranchAddress("run",&runNumber);
  cout<<"Opened the file and set the Branches"<<endl; 
  
  ///Get the total number of entries within the tree
  Int_t nument=ceventTree->GetEntries();

  vector <Double_t> UnixTime[MCH];
  vector <Double_t> V2n[MCH];
  vector <Double_t> CorScore[MCH];
  vector <Double_t> CorTime[MCH];
  vector <Double_t> V2nAvg;
  vector <Double_t> CorScoreAvg;

  ceventTree->GetEntry(0);

  TString PrintFile="ARA";
  PrintFile+=StationId;
  PrintFile+="_Run";
  PrintFile+=runNumber;
  PrintFile+=".txt";

  const char* filename=PrintFile.Data();

  ofstream aout(filename);
  
  TH2D * h2=new TH2D("","",200,0,60000,200,0,12000);
  
  TGraph *grCor[MCH];
  TGraph *gr[MCH];

  const Int_t NumCutCh=0;
  Double_t CutCh[NumCutCh]={};
  Double_t FirstUnixTime=0;
  
  ///Start looping over the events
  for(Int_t ievt=5; ievt<nument; ievt++){
  
    ceventTree->GetEntry(ievt);
    cout<<ievt<<" "<<rawAtriEvPtr->eventNumber<<" "<<rawAtriEvPtr->isCalpulserEvent()<<endl;
    ///////Right now only look at RF events which are not calpulser and software triggers
    //if(rawAtriEvPtr->isCalpulserEvent()==false && rawAtriEvPtr->timeStamp>1000 && rawAtriEvPtr->isSoftwareTrigger()==false){
      ///Initialise the Useful event pointer to load the ARA event waveforms with the right calibration
      //UsefulAtriStationEvent *realAtriEvPtr=new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
      
      // Double_t EventNum=rawAtriEvPtr->eventNumber;
      // Double_t UTime=rawAtriEvPtr->unixTime;
      // Double_t UTimeUs=rawAtriEvPtr->unixTimeUs;
      // Double_t TStamp=rawAtriEvPtr->timeStamp;
      // if(ievt==5){
      // 	FirstUnixTime=UTime;
      // 	cout<<"UTime "<<setprecision(10)<<UTime<<endl;
      // }
      // // vector <Double_t> yar[MCH];   
      // // vector <Double_t> xar[MCH];   

      // Double_t V2nSum=0;
      // Double_t CorScoreSum=0;
      
      // ///For each event loop over the 16 channels
      // for(Int_t ich=0; ich<MCH; ich++){
      // 	Bool_t SkipCh=false;
      // 	for(Int_t icut=0; icut<NumCutCh; icut++){
      // 	  if(ich==CutCh[icut]){
      // 	    SkipCh=true;
      // 	  }
      // 	}
      // 	if(SkipCh==false){
      // 	  ///Get the Waveform from the data file for each channel
      // 	  TGraph *grdum=realAtriEvPtr->getGraphFromRFChan(ich);

      // 	  ///Interpolate the waveforms to ensure equal spacing between the samples
      // 	  TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.6);
      // 	  delete grdum;
	
      // 	  ///Get the SNR of the waveform of each channel
      // 	  //cout<<"SNR of channel "<<ich<<" is "<<FFTtools::getWaveformSNR(gr[ich])<<endl;
	
      // 	  ///Get the total no. of points in the waveform
      // 	  Int_t grpnts=gr->GetN();	

      // 	  Double_t V2=0;
      // 	  for(Int_t i=0; i<grpnts; i++){
      // 	    Double_t x,y;
      // 	    gr->GetPoint(i,x,y);
      // 	    // xar[ich].push_back(x);
      // 	    // yar[ich].push_back(y*y);
      // 	    V2+=y*y;
      // 	  }
      // 	  V2nSum+=V2/(double)grpnts;
      // 	  UnixTime[ich].push_back(UTime-FirstUnixTime);
      // 	  V2n[ich].push_back(V2/(double)grpnts);
	  
      // 	  if(ich<8){
      // 	    grCor[ich]=FFTtools::getCorrelationGraph(gCPtemp[0][ich],gr);
      // 	  }else{
      // 	    grCor[ich]=FFTtools::getCorrelationGraph(gCPtemp[1][ich],gr);
      // 	  }
      // 	  Int_t nBins = grCor[ich]->GetN();
      // 	  Double_t *yVals = grCor[ich]->GetY();
      // 	  Double_t *xVals = grCor[ich]->GetX();
	  
      // 	  Double_t MaxCorScore=TMath::MaxElement(nBins,yVals);
      // 	  Int_t CorTimeBin=TMath::LocMax(nBins,yVals);
      // 	  Double_t CorTimeVal=xVals[CorTimeBin];

      // 	  CorScoreSum+=MaxCorScore;
      // 	  CorScore[ich].push_back(MaxCorScore);
      // 	  CorTime[ich].push_back(CorTimeVal);
	  
      // 	  //cout<<ich<<" "<<CorScore[ich+MCH*itmp]<<" "<<CorTime[ich+MCH*itmp]<<endl;
	  
      // 	  delete []xVals;
      // 	  delete []yVals;
	
      // 	  delete gr;
      // 	}
      // }////channel loop
      
      // V2nAvg.push_back(V2nSum/(MCH-NumCutCh));
      // CorScoreAvg.push_back(CorScoreSum/(MCH-NumCutCh));
      // h2->Fill(V2nSum/(MCH-NumCutCh),CorScoreSum/(MCH-NumCutCh));
      
      // //cout<<ievt<<" "<<V2nSum/(MCH-NumCutCh)<<" "<<CorScoreSum/(MCH-NumCutCh)<<endl;
      // if(V2nSum/(MCH-NumCutCh)>2000 && CorScoreSum/(MCH-NumCutCh)>1000){
      // 	aout<<runNumber<<" "<<ievt<<endl;
      // 	//cout<<"It is a SPICE event! "<<runNumber<<" "<<ievt<<endl;
      //  }
    //delete realAtriEvPtr;      
      aout<<runNumber<<" "<<ievt<<" "<<rawAtriEvPtr->eventNumber<<endl;
      // }
  
  }
  
  // TMultiGraph * mgV2=new TMultiGraph();
  // TMultiGraph * mgCorScore=new TMultiGraph(); 
  // TGraph *grV2[MCH];
  // TGraph *grCorScore[MCH];
  
  // for(Int_t ich=0; ich<MCH; ich++){
  //   Bool_t SkipCh=false;
  //   for(Int_t icut=0; icut<NumCutCh; icut++){
  //     if(ich==CutCh[icut]){
  // 	SkipCh=true;
  //     }
  //   }
  //   if(SkipCh==false){

  //     grV2[ich]=new TGraph(V2n[ich].size(),UnixTime[ich].data(),V2n[ich].data());
  //     grCorScore[ich]=new TGraph(CorScore[ich].size(),UnixTime[ich].data(),CorScore[ich].data());
    
  //     TString title="Ch. ";
  //     title+=ich;
  //     grV2[ich]->SetTitle(title);
  //     grCorScore[ich]->SetTitle(title);
    
  //     grV2[ich]->SetMarkerStyle(20);
  //     grCorScore[ich]->SetMarkerStyle(20);
 
  //     Int_t Color=0;
  //     if(ich<8){
  // 	Color=1+ich;
  //     }
    
  //     if(ich==9){
  // 	Color=13;
  //     }
  //     if(ich==10){
  // 	Color=28;
  //     }
  //     if(ich==12){
  // 	Color=20;
  //     }
  //     if(ich==14){
  // 	Color=30;
  //     }
  //     grV2[ich]->SetMarkerColor(Color);
  //     grCorScore[ich]->SetMarkerColor(Color);
    
  //     mgV2->Add(grV2[ich]);
  //     mgCorScore->Add(grCorScore[ich]);
    
  //   }
  // }

  // TGraph *grV2nAvg=new TGraph(V2nAvg.size(),UnixTime[0].data(),V2nAvg.data());
  // TGraph *grCorScoreAvg=new TGraph(CorScoreAvg.size(),UnixTime[0].data(),CorScoreAvg.data());
  // grV2nAvg->SetMarkerStyle(20);
  // grCorScoreAvg->SetMarkerStyle(20);   

  // grV2nAvg->SetMarkerColor(kRed);
  // grCorScoreAvg->SetMarkerColor(kBlue);
  
  // TCanvas *c1 = new TCanvas("c1","c1",1800,1800);
  // c1->cd();
  // c1->SetGridx();
  // c1->SetGridy();
  // mgV2->GetXaxis()->SetLabelSize(0.05);
  // mgV2->GetYaxis()->SetLabelSize(0.05);
  // mgV2->GetXaxis()->SetTitleSize(0.05);
  // mgV2->GetYaxis()->SetTitleSize(0.05);

  // mgV2->Draw("AP");
  // mgV2->SetTitle(";Unixtime (s);Avg. Sample Power in Channel;");
  // c1->BuildLegend();
  // c1->SaveAs("AvgPower_AllCh_run12576.png");
  // c1->SaveAs("AvgPower_AllCh_run12576.pdf");
    
  // TCanvas *c2 = new TCanvas("c2","c2",1800,1800);
  // c2->cd();
  // c2->SetGridx();
  // c2->SetGridy();
  // mgCorScore->GetXaxis()->SetLabelSize(0.05);
  // mgCorScore->GetYaxis()->SetLabelSize(0.05);
  // mgCorScore->GetXaxis()->SetTitleSize(0.05);
  // mgCorScore->GetYaxis()->SetTitleSize(0.05);
  // mgCorScore->Draw("AP");
  // mgCorScore->SetTitle(";Unixtime (s);Peak Correlation Score;");
  // c2->BuildLegend();
  // c2->SaveAs("PeakCorScore_AllCh_run12576.png");
  // c2->SaveAs("PeakCorScore_AllCh_run12576.pdf");
  
  // TCanvas *c3 = new TCanvas("c3","c3",1800,1800);
  // c3->Divide(1,2);
  // c3->cd(1);
  // c3->cd(1)->SetGridx();
  // c3->cd(1)->SetGridy();
  // grV2nAvg->GetXaxis()->SetLabelSize(0.05);
  // grV2nAvg->GetYaxis()->SetLabelSize(0.05);
  // grV2nAvg->GetXaxis()->SetTitleSize(0.05);
  // grV2nAvg->GetYaxis()->SetTitleSize(0.05);
  // grV2nAvg->Draw("AP");
  // grV2nAvg->SetTitle(";Unixtime (s);Avg. Sample Power Over All Channels;");
  
  // c3->cd(2);
  // c3->cd(2)->SetGridx();
  // c3->cd(2)->SetGridy();
  // grCorScoreAvg->GetXaxis()->SetLabelSize(0.05);
  // grCorScoreAvg->GetYaxis()->SetLabelSize(0.05);
  // grCorScoreAvg->GetXaxis()->SetTitleSize(0.05);
  // grCorScoreAvg->GetYaxis()->SetTitleSize(0.05);
  // grCorScoreAvg->Draw("AP");
  // grCorScoreAvg->SetTitle(";Unixtime (s);Avg. Peak Correlation Score Over All Channels;");
  // c3->SaveAs("AvgPowerNCorScore_AllCh_run12576.png");
  // c3->SaveAs("AvgPowerNCorScore_AllCh_run12576.pdf");
 
  // TCanvas *c4 = new TCanvas("c4","c4",1800,1800);
  // c4->cd();
  // c4->SetGridx();
  // c4->SetGridy();
  // h2->GetXaxis()->SetLabelSize(0.05);
  // h2->GetYaxis()->SetLabelSize(0.05);
  // h2->GetXaxis()->SetTitleSize(0.05);
  // h2->GetYaxis()->SetTitleSize(0.05);
  // h2->Draw("colz");
  // h2->SetTitle(";Avg. Sample Power Over All Channels;Avg. Peak Correlation Score Over All Channels;"); 
  // c4->SaveAs("AvgPowerNCorScore_AllCh_colz_run12576.png");
  //c4->SaveAs("AvgPowerNCorScore_AllCh_colz_run12576.pdf");

}
