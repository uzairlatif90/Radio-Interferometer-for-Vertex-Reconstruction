const Int_t MCH=16;

#include "FFTtools.h"
#include "../Interferometer.cc"

Double_t TrueX=-456.721;
Double_t TrueY=-2353;

TGraph *gCPtemp[2][16];

void ReadCPTemp(){  

  {
    TString filename="";
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
    TString filename="";
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

double GetNoiseRMS(TGraph *gr){

  int nDiv=10;
  int nPoints = gr->GetN();
  double *yVals = gr->GetY();
  int PPD=floor((double)nPoints/nDiv);/// divide the rest of the WF into 10 parts
    
  double mean=0.;
  double meanSq=0.;
  ///find the part with the lowest RMS
  double newRMS=10000000;
  double rms=0;
  int noidiv=0;
  for(int nd=0;nd<nDiv;nd++){
    mean=0.;
    meanSq=0.;
    for(int i=nd*PPD;i<PPD+(nd*PPD);i++){
      mean+=yVals[i];
      meanSq+=yVals[i]*yVals[i];
    }
    mean=mean/PPD;
    meanSq=meanSq/PPD;
    rms=sqrt(meanSq-mean*mean);
    if(rms<newRMS){
      newRMS=rms;
      noidiv=nd;
    }
  }
  
  return newRMS;
}

void PeakFinder(TGraph *grPwrEnvOriginal, TGraph *grPeakPoint){

  TGraph *grPwrEnvSmooth=FFTtools::getBoxCar(grPwrEnvOriginal,7); 
 	
  vector <double> Peaks[3];
  vector <double> Troughs[3];
	
  int TotalgrSamples=grPwrEnvSmooth->GetN();
  int DeltaSamples=8;
  double yi=0;
  for(int isample=0; isample<TotalgrSamples; isample++){
    double xp,yp;
    grPwrEnvSmooth->GetPoint(isample,xp,yp);

    double RightSideCount=0;
    for(int isampleDelta=isample; isampleDelta<=DeltaSamples+isample; isampleDelta++){
      if(isampleDelta<TotalgrSamples){
	grPwrEnvSmooth->GetPoint(isampleDelta,xp,yp);
	if(isampleDelta==isample){
	  yi=yp;
	}
	if(yp<yi){
	  RightSideCount++;
	}
      }else{
	isampleDelta=DeltaSamples+isample+1;
      }
    }

    double LeftSideCount=0;
    for(int isampleDelta=isample; isampleDelta>=isample-DeltaSamples; isampleDelta--){
      if(isampleDelta>=0){
	grPwrEnvSmooth->GetPoint(isampleDelta,xp,yp);
	if(isampleDelta==isample){
	  yi=yp;
	}
	if(yp<yi){
	  LeftSideCount++;
	}
      }else{
	isampleDelta=isample-DeltaSamples-1;
      }
    }

    grPwrEnvSmooth->GetPoint(isample,xp,yp);
    if(LeftSideCount==DeltaSamples && RightSideCount==DeltaSamples){
      //cout<<"we have a peak "<<xp<<" "<<yp<<endl;
      Peaks[0].push_back(xp);
      Peaks[1].push_back(yp);
      Peaks[2].push_back(isample);
    }

    if(LeftSideCount+RightSideCount<=DeltaSamples*2 && isample<=DeltaSamples*2+1 && RightSideCount>0 && LeftSideCount>0){
      //cout<<"we have a peak "<<xp<<" "<<yp<<endl;
      Peaks[0].push_back(xp);
      Peaks[1].push_back(yp);
      Peaks[2].push_back(isample);
    }
    if(RightSideCount==DeltaSamples  && isample<=DeltaSamples && LeftSideCount==0){
      //cout<<"we have a peak "<<xp<<" "<<yp<<endl;
      Peaks[0].push_back(xp);
      Peaks[1].push_back(yp);
      Peaks[2].push_back(isample);
    }

    if(LeftSideCount+RightSideCount<=DeltaSamples*2 && isample>=TotalgrSamples-(DeltaSamples*2+1) && RightSideCount>0 && LeftSideCount>0){
      //cout<<"we have a peak "<<xp<<" "<<yp<<endl;
      Peaks[0].push_back(xp);
      Peaks[1].push_back(yp);
      Peaks[2].push_back(isample);
    }
    if(LeftSideCount==DeltaSamples  && isample>=TotalgrSamples-DeltaSamples && RightSideCount==0){
      //cout<<"we have a peak "<<xp<<" "<<yp<<endl;
      Peaks[0].push_back(xp);
      Peaks[1].push_back(yp);
      Peaks[2].push_back(isample);
    }
	  
  }

  //int DummyBin=TMath::LocMax(TotalgrSamples,grPwrEnvSmooth->GetY());
  // if(DummyBin<DeltaSamples || DummyBin>TotalgrSamples-DeltaSamples){
  //   double xp,yp;
  //   grPwrEnvSmooth->GetPoint(DummyBin,xp,yp);
  //   Peaks[0].push_back(xp);
  //   Peaks[1].push_back(yp);
  //   Peaks[2].push_back(DummyBin);
  // }
  
  double PowerPeakAmp[4]={0,0,0,0};
  double PowerPeakTime[4]={0,0,0,0};
  double PowerPeakBin[4]={0,0,0,0};
	
  int DummyBin=TMath::LocMax(Peaks[1].size(),Peaks[1].data());
  PowerPeakAmp[0]=Peaks[1][DummyBin];
  PowerPeakTime[0]=Peaks[0][DummyBin];
  PowerPeakBin[0]=Peaks[2][DummyBin];
  Peaks[1][DummyBin]=0;
	
  DummyBin=TMath::LocMax(Peaks[1].size(),Peaks[1].data());
  PowerPeakAmp[1]=Peaks[1][DummyBin];
  PowerPeakTime[1]=Peaks[0][DummyBin];
  PowerPeakBin[1]=Peaks[2][DummyBin];

  for(int ipeak=0;ipeak<2;ipeak++){
    vector <double> RefinePeak[3];
    for(int isample=PowerPeakBin[ipeak]-DeltaSamples; isample<=PowerPeakBin[ipeak]+DeltaSamples; isample++){
      double xp,yp;
      grPwrEnvOriginal->GetPoint(isample,xp,yp);
      RefinePeak[0].push_back(xp);
      RefinePeak[1].push_back(yp);
      RefinePeak[2].push_back(isample);
    }
      
    DummyBin=TMath::LocMax(RefinePeak[1].size(),RefinePeak[1].data());
    PowerPeakAmp[ipeak]=RefinePeak[1][DummyBin];
    PowerPeakTime[ipeak]=RefinePeak[0][DummyBin];
    PowerPeakBin[ipeak]=RefinePeak[2][DummyBin];
    RefinePeak[0].clear();
    RefinePeak[1].clear();
    RefinePeak[2].clear();
  }  

  double NoiseRMS=GetNoiseRMS(grPwrEnvOriginal);
  int IgnorePeak[2]={1,1};
  double LargePeak=PowerPeakAmp[1];
  double SmallPeak=PowerPeakAmp[0];
  if(SmallPeak>LargePeak){
    swap(SmallPeak,LargePeak);
  }

  if(SmallPeak>=LargePeak*0.25 && fabs(PowerPeakTime[1]-PowerPeakTime[0])>40 ){
    //cout<<"We have two peaks "<<endl;

    if(PowerPeakTime[0]>PowerPeakTime[1]){
      swap(PowerPeakTime[0],PowerPeakTime[1]);
      swap(PowerPeakAmp[0],PowerPeakAmp[1]);
      swap(PowerPeakBin[0],PowerPeakBin[1]);
    }

    for(int ipeak=0;ipeak<2;ipeak++){
      vector <double> RefinePeak[3];
   
      for(int isample=PowerPeakBin[ipeak]-20; isample<=PowerPeakBin[ipeak]+20; isample++){
	if(isample>-1){
	  double xp,yp;
	  grPwrEnvOriginal->GetPoint(isample,xp,yp);
	  if(isample<PowerPeakBin[ipeak]-5 || isample>PowerPeakBin[ipeak]+5){
	    RefinePeak[0].push_back(xp);
	    RefinePeak[1].push_back(yp);
	    RefinePeak[2].push_back(isample);
	  }
	}
      }
      
      DummyBin=TMath::LocMax(RefinePeak[1].size(),RefinePeak[1].data());
      PowerPeakAmp[ipeak+2]=RefinePeak[1][DummyBin];
      PowerPeakTime[ipeak+2]=RefinePeak[0][DummyBin];
      PowerPeakBin[ipeak+2]=RefinePeak[2][DummyBin];
      RefinePeak[0].clear();
      RefinePeak[1].clear();
      RefinePeak[2].clear();   
    
      double LargePeak=PowerPeakAmp[ipeak+2];
      double SmallPeak=PowerPeakAmp[ipeak];
      if(SmallPeak>LargePeak){
	swap(SmallPeak,LargePeak);
      }
      if(fabs(PowerPeakTime[ipeak+2]-PowerPeakTime[ipeak])>25*0.6 && SmallPeak>=LargePeak*0.10){
	PowerPeakAmp[ipeak]=(PowerPeakAmp[ipeak]+PowerPeakAmp[ipeak+2])/2;
	PowerPeakTime[ipeak]=(PowerPeakTime[ipeak]*PowerPeakAmp[ipeak]+PowerPeakTime[ipeak+2]*PowerPeakAmp[ipeak+2])/(PowerPeakAmp[ipeak]+PowerPeakAmp[ipeak+2]);
      }
    }
    
    if(PowerPeakBin[1]<2*DeltaSamples){
      IgnorePeak[1]=0;
    }
    if(PowerPeakBin[0]<2*DeltaSamples){
      IgnorePeak[0]=0;
    }
    
    // int StartBin=PowerPeakBin[0];
    // int EndBin=PowerPeakBin[1];
    // if(StartBin>EndBin){
    //   StartBin=PowerPeakBin[1];
    //   EndBin=PowerPeakBin[0];
    // }
    // for(int isample=StartBin; isample<EndBin; isample++){
    //   double xp,yp;
    //   grPwrEnvOriginal->GetPoint(isample,xp,yp);
    //   Troughs[0].push_back(xp);
    //   Troughs[1].push_back(yp);
    //   Troughs[2].push_back(isample);
    // }
	  
    // DummyBin=TMath::LocMin(Troughs[1].size(),Troughs[1].data());
    // double TroughAmp=Troughs[1][DummyBin];
    // double TroughTime=Troughs[0][DummyBin];
    // double TroughBin=Troughs[2][DummyBin];

    //cout<<PowerPeakTime[0]<<" "<<PowerPeakAmp[0]<<" "<<PowerPeakTime[1]<<" "<<PowerPeakAmp[1]<<" "<<TroughTime<<" "<<TroughAmp<<endl;
  
    grPeakPoint->SetPoint(0,PowerPeakTime[0],PowerPeakAmp[0]);
    //grPeakPoint->SetPoint(1,TroughTime,TroughAmp);
    grPeakPoint->SetPoint(1,PowerPeakTime[1],PowerPeakAmp[1]);
    grPeakPoint->SetPoint(2,PowerPeakTime[0],PowerPeakAmp[0]/NoiseRMS);
    grPeakPoint->SetPoint(3,PowerPeakTime[1],PowerPeakAmp[1]/NoiseRMS);
    grPeakPoint->SetPoint(4,IgnorePeak[0],0);
    grPeakPoint->SetPoint(5,IgnorePeak[1],0);
  }else{
    //cout<<"We have one peak "<<endl;

    for(int ipeak=0;ipeak<1;ipeak++){
      vector <double> RefinePeak[3];
   
      for(int isample=PowerPeakBin[ipeak]-20; isample<=PowerPeakBin[ipeak]+20; isample++){
	if(isample>-1){
	  double xp,yp;
	  grPwrEnvOriginal->GetPoint(isample,xp,yp);
	  if(isample<PowerPeakBin[ipeak]-5 || isample>PowerPeakBin[ipeak]+5){
	    RefinePeak[0].push_back(xp);
	    RefinePeak[1].push_back(yp);
	    RefinePeak[2].push_back(isample);
	  }
	}
      }
      
      DummyBin=TMath::LocMax(RefinePeak[1].size(),RefinePeak[1].data());
      PowerPeakAmp[ipeak+2]=RefinePeak[1][DummyBin];
      PowerPeakTime[ipeak+2]=RefinePeak[0][DummyBin];
      PowerPeakBin[ipeak+2]=RefinePeak[2][DummyBin];
      RefinePeak[0].clear();
      RefinePeak[1].clear();
      RefinePeak[2].clear(); 
    
      double LargePeak=PowerPeakAmp[ipeak+2];
      double SmallPeak=PowerPeakAmp[ipeak];
      if(SmallPeak>LargePeak){
	swap(SmallPeak,LargePeak);
      }
      if(fabs(PowerPeakTime[ipeak+2]-PowerPeakTime[ipeak])>25*0.6 && SmallPeak>=LargePeak*0.10){
	PowerPeakAmp[ipeak]=(PowerPeakAmp[ipeak]+PowerPeakAmp[ipeak+2])/2;
	PowerPeakTime[ipeak]=(PowerPeakTime[ipeak]*PowerPeakAmp[ipeak]+PowerPeakTime[ipeak+2]*PowerPeakAmp[ipeak+2])/(PowerPeakAmp[ipeak]+PowerPeakAmp[ipeak+2]);
      }
    }
    
    if(PowerPeakBin[0]<2*DeltaSamples){
      IgnorePeak[0]=0;
    }
    
    //cout<<PowerPeakTime[0]<<" "<<PowerPeakAmp[0]<<" "<<PowerPeakTime[1]<<" "<<PowerPeakAmp[1]<<endl;
    grPeakPoint->SetPoint(0,PowerPeakTime[0],PowerPeakAmp[0]);
    grPeakPoint->SetPoint(1,PowerPeakTime[0],PowerPeakAmp[0]/NoiseRMS);
    grPeakPoint->SetPoint(2,IgnorePeak[0],0);
  }
      
  Peaks[0].clear();
  Peaks[1].clear();
  Peaks[2].clear();

  Troughs[0].clear();
  Troughs[1].clear();
  Troughs[2].clear();
  
}

void GetSPICEWaveform(){

  ReadCPTemp();
  Int_t StationId=2;
  DeclareAntennaConfigARA(StationId);

  Double_t Event=3900;
  
  ////initialise the event pointer to take data from the event tree
  RawAtriStationEvent *rawAtriEvPtr=0;

  ////initialise the AraEventCalibrator class  
  AraEventCalibrator *calib = AraEventCalibrator::Instance();
  
  ///Open the ARA root file that contains the data
  TFile *newfile=new TFile("../../../ARA_Stuff/TestInterSPICEdata/run_012576/event012576.root");
   
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
  
  double ChHitTime[2][TotalAntennasRx];
  int IgnoreCh[2][TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];	
  double ChAmp[2][TotalAntennasRx];
  
  TGraph *grCor[MCH];
  TGraph *gr[MCH];
  TGraph *grPwrEnv[MCH];
  TGraph *grWF[MCH];
  
  const Int_t NumCutCh=4;
  Double_t CutCh[NumCutCh]={8,15,11,13};
  Double_t FirstUnixTime=0;
  
  ///Start looping over the events
  for(Int_t ievt=Event; ievt<Event+100; ievt++){
  
    ceventTree->GetEntry(ievt);
    cout<<ievt<<" "<<rawAtriEvPtr->eventNumber<<" "<<rawAtriEvPtr->unixTime-1545775998<<endl;
    /////Right now only look at RF events which are not calpulser and software triggers
    if(rawAtriEvPtr->isCalpulserEvent()==false && rawAtriEvPtr->timeStamp>1000 && rawAtriEvPtr->isSoftwareTrigger()==false){
      ///Initialise the Useful event pointer to load the ARA event waveforms with the right calibration
      UsefulAtriStationEvent *realAtriEvPtr=new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
      
      Double_t EventNum=rawAtriEvPtr->eventNumber;
      Double_t UTime=rawAtriEvPtr->unixTime;
      Double_t UTimeUs=rawAtriEvPtr->unixTimeUs;
      Double_t TStamp=rawAtriEvPtr->timeStamp;
      if(ievt==5){
	FirstUnixTime=UTime;
	cout<<"UTime "<<setprecision(10)<<UTime<<endl;
      }
      // vector <Double_t> yar[MCH];   
      // vector <Double_t> xar[MCH];   

      Double_t V2nSum=0;
      Double_t CorScoreSum=0;
      
      ///For each event loop over the 16 channels
      for(Int_t ich=0; ich<MCH; ich++){
	Bool_t SkipCh=false;
	for(Int_t icut=0; icut<NumCutCh; icut++){
	  if(ich==CutCh[icut]){
	    SkipCh=true;
	  }
	}
	if(SkipCh==false){
	  ///Get the Waveform from the data file for each channel
	  TGraph *grdum=realAtriEvPtr->getGraphFromRFChan(ich);

	  ///Interpolate the waveforms to ensure equal spacing between the samples
	  TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.6);
	  delete grdum;
	
	  ///Get the SNR of the waveform of each channel
	  //cout<<"SNR of channel "<<ich<<" is "<<FFTtools::getWaveformSNR(gr[ich])<<endl;
	
	  ///Get the total no. of points in the waveform
	  Int_t grpnts=gr->GetN();	

	  Double_t V2=0;
	  for(Int_t i=0; i<grpnts; i++){
	    Double_t x,y;
	    gr->GetPoint(i,x,y);
	    // xar[ich].push_back(x);
	    // yar[ich].push_back(y*y);
	    V2+=y*y;
	  }
	  V2nSum+=V2/(double)grpnts;
	  UnixTime[ich].push_back(UTime-FirstUnixTime);
	  V2n[ich].push_back(V2/(double)grpnts);
	  
	  if(ich<8){
	    grCor[ich]=FFTtools::getCorrelationGraph(gCPtemp[0][ich],gr);
	  }else{
	    grCor[ich]=FFTtools::getCorrelationGraph(gCPtemp[1][ich],gr);
	  }
	  Int_t nBins = grCor[ich]->GetN();
	  Double_t *yVals = grCor[ich]->GetY();
	  Double_t *xVals = grCor[ich]->GetX();
	  
	  Double_t MaxCorScore=TMath::MaxElement(nBins,yVals);
	  Int_t CorTimeBin=TMath::LocMax(nBins,yVals);
	  Double_t CorTimeVal=xVals[CorTimeBin];

	  CorScoreSum+=MaxCorScore;
	  CorScore[ich].push_back(MaxCorScore);
	  CorTime[ich].push_back(CorTimeVal);
	  
	  //cout<<ich<<" "<<CorScore[ich+MCH*itmp]<<" "<<CorTime[ich+MCH*itmp]<<endl;
	  
	  delete []xVals;
	  delete []yVals;
	
	  delete gr;
	}
      }////channel loop
      
      V2nAvg.push_back(V2nSum/(MCH-NumCutCh));
      CorScoreAvg.push_back(CorScoreSum/(MCH-NumCutCh));
      
      //cout<<ievt<<" "<<V2nSum/(MCH-NumCutCh)<<" "<<CorScoreSum/(MCH-NumCutCh)<<endl;
      if(V2nSum/(MCH-NumCutCh)>2000 && CorScoreSum/(MCH-NumCutCh)>1000 && UTime-1545775998<1000){
	cout<<runNumber<<" "<<ievt<<endl;
	cout<<"It is a SPICE event! "<<runNumber<<" "<<ievt<<endl;

	double dtDR_Avg=0;
	int IgnorePeakCount=0;
	for(int ich=0; ich<MCH; ich++){
	  Bool_t SkipCh=false;
	  for(int icut=0; icut<NumCutCh; icut++){
	    if(ich==CutCh[icut]){
	      SkipCh=true;
	    }
	  }
	  if(SkipCh==false){
	    ///Get the Waveform from the data file for each channel
	    TGraph *grdum=realAtriEvPtr->getGraphFromRFChan(ich);

	    ///Interpolate the waveforms to ensure equal spacing between the samples
	    TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.6);
	    delete grdum;

	    grPwrEnv[ich]=FFTtools::getSimplePowerEnvelopeGraph(gr);
	    TGraph *grPeakPoint=new TGraph();
	    //grWF[ich]=new TGraph();
	    grWF[ich]=FFTtools::getSimplePowerEnvelopeGraph(gr);
	   
	    PeakFinder(grPwrEnv[ich], grPeakPoint);
	  
	    int NumOfPulses=grPeakPoint->GetN();
	    double Dtime=0,Rtime=0;
	    double DAmp=0,RAmp=0;
	    double DSNR=0,RSNR=0;
	    int IgnorePeakD=0,IgnorePeakR=0;
	    if(NumOfPulses==6){
	      double xpeak,ypeak;
	      grPeakPoint->GetPoint(0,xpeak,ypeak);
	      Dtime=xpeak;
	      DAmp=ypeak;
	      grPeakPoint->GetPoint(2,xpeak,ypeak);
	      DSNR=ypeak;
	      grPeakPoint->GetPoint(4,xpeak,ypeak);
	      IgnorePeakD=xpeak;
	   
	      grPeakPoint->GetPoint(1,xpeak,ypeak);
	      Rtime=xpeak;
	      RAmp=ypeak;
	      grPeakPoint->GetPoint(3,xpeak,ypeak);
	      RSNR=ypeak;
	      grPeakPoint->GetPoint(5,xpeak,ypeak);
	      IgnorePeakR=xpeak;
	    
	      if(Dtime>Rtime){
		swap(Dtime,Rtime);
		swap(DAmp,RAmp);
		swap(DSNR,RSNR);
		swap(IgnorePeakD,IgnorePeakR);
	      }
	      dtDR_Avg+=fabs(Dtime-Rtime);

	      if(IgnorePeakD==0){
		Dtime=0;
		DSNR=0;
		DAmp=0;
		IgnorePeakCount++;
	      }

	      if(IgnorePeakR==0){
		Rtime=0;
		RSNR=0;
		RAmp=0;
		IgnorePeakCount++;
	      }
	      
	      ChHitTime[0][ich]=Dtime;
	      ChHitTime[1][ich]=Rtime;
	      ChSNR[0][ich]=DSNR;
	      ChSNR[1][ich]=RSNR;
	      ChAmp[0][ich]=DAmp;
	      ChAmp[1][ich]=RAmp;
	    }else{
	      double xpeak,ypeak;
	      grPeakPoint->GetPoint(0,xpeak,ypeak);
	      Dtime=xpeak;
	      DAmp=ypeak;
	      grPeakPoint->GetPoint(1,xpeak,ypeak);
	      DSNR=ypeak;
	      grPeakPoint->GetPoint(2,xpeak,ypeak);
	      IgnorePeakD=xpeak;
	    
	      if(IgnorePeakD==0){
		Dtime=0;
		DSNR=0;
		DAmp=0;
		IgnorePeakCount++;
	      }
	      ChHitTime[0][ich]=Dtime;
	      ChHitTime[1][ich]=0;
	      ChSNR[0][ich]=DSNR;
	      ChSNR[1][ich]=0;
	      ChAmp[0][ich]=DAmp;
	      ChAmp[1][ich]=0;
	    }
	  
	    delete gr;
	    delete grPeakPoint;
	  }
	}////channel loop
      break;	
      }
      delete realAtriEvPtr;      

    }
  }

  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(4,4);
  for(int ich=0; ich<MCH; ich++){
    Bool_t SkipCh=false;
    for(int icut=0; icut<NumCutCh; icut++){
      if(ich==CutCh[icut]){
  	SkipCh=true;
      }
    }
    if(SkipCh==false){
      c1->cd(ich+1);
      grWF[ich]->Draw("AL");
    }
  }
  
  
}
