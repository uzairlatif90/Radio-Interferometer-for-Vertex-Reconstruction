const int MCH=16;

#include "FFTtools.h"
#include "../../../RET_int/Interferometer/Interferometer.cc"
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

gsl_interp_accel *spline_acc;
gsl_spline *spline_steffen;

void LoadDepthFile(){
  
  int StartTime=1545822000-13*60*60;
  ////Open the file
  std::ifstream ain("../Depth26Dec.txt");
  int n1=0;////variable for counting total number of data points
  std::string line;
  int dummy[2]={0,0};////temporary variable for storing data values from the file  

  vector<double> Time[3];
  vector<double> Depth;
  //Check if file is open and store data
  if(ain.is_open()){
    while (true){
      ain>>dummy[0]>>dummy[1];
      if( ain.eof()) break;
      double hrs=(dummy[0]-dummy[0]%100)/100 -11;
      double min=dummy[0]%100;
      Time[0].push_back(hrs);
      Time[1].push_back(min);
      Time[2].push_back(StartTime+min*60+hrs*60*60);
      Depth.push_back(dummy[1]);
      //cout<<setprecision(10)<<n1<<" "<<hrs<<" "<<min<<" "<<StartTime+min*60+hrs*60*60<<" "<<dummy[1]<<endl;
      n1++;
    }
  }
  ain.close();

  spline_acc = gsl_interp_accel_alloc();
  spline_steffen = gsl_spline_alloc(gsl_interp_steffen, Time[0].size());
  gsl_spline_init(spline_steffen, Time[2].data(), Depth.data(), Time[0].size());

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

void ReconstructSPICEevents(char const *InputFileName, int Run, int Event){

  TString OutputFileName="./output/";
  OutputFileName+="Run";
  OutputFileName+=Run;
  OutputFileName+="Event";
  OutputFileName+=Event;
  OutputFileName+=".root";

  int StationId=2;
  const int NumCutCh=4;
  double CutCh[NumCutCh]={8,15,11,13};
  double TrueX=-456.721;
  double TrueY=-2353;
  double ExpectedPositionUncertainty=5;//in m
  
  DeclareAntennaConfigARA(StationId);
  LoadDepthFile();
  
  ////initialise the event pointer to take data from the event tree
  RawAtriStationEvent *rawAtriEvPtr=0;

  ////initialise the AraEventCalibrator class  
  AraEventCalibrator *calib = AraEventCalibrator::Instance();
  
  ///Open the ARA root file that contains the data
  // TFile *newfile=new TFile("../run_012576/event012576.root");
  // TFile *ResultFile=new TFile("Reco12576_results.root","RECREATE");

  TFile *InputFile=new TFile(InputFileName);
  
  ///Open the Tree in the file that contains data from all of the events
  TTree *ceventTree = (TTree*)InputFile->Get("eventTree");

  int runNumber;
  ///Set the Branch addresses
  ceventTree->SetBranchAddress("event",&rawAtriEvPtr);
  ceventTree->SetBranchAddress("run",&runNumber);
  cout<<"Opened the file and set the Branches"<<endl; 
  
  ///Get the total number of entries within the tree
  int nument=ceventTree->GetEntries();

  vector <double> UnixTime[MCH];
  vector <double> V2n[MCH];
  vector <double> PwrSNR[2][MCH];
  
  vector <double> dtDR_Ch[MCH];
  vector <double> DRAmpRatio_Ch[MCH];
  vector <double> SPICEReco[3];
  vector <double> UnixTimeSelected;
  vector <double> UnixTimeReco;
  vector <double> DurationReco;

  vector <double> V2nAvg;
  vector <double> CorScoreAvg;

  TGraph *gr[MCH];
  TGraph *grPwrEnv[MCH];

  double eventNum;
  double unixTime;
  double firstUnixTime;
  double timeStamp;
  int runNum;
  double Duration;
  double FinalMinValue;
  int Iterations;
  double FinalTxCor_XYZ[3];
  double FinalTxCor_ThPhR[3];
  double InitialTxCor_XYZ[3];
  double InitialTxCor_ThPhR[3];
  double dX,dY,dZ;
  double dTh,dPh,dR;
  
  TFile *OutputFile=new TFile(OutputFileName,"RECREATE");  

  TTree *RecoTree = new TTree("RecoTree","Reco info about Event");
  RecoTree->Branch("eventNum",&eventNum,"eventNum/D");
  RecoTree->Branch("unixTime",&unixTime,"unixTime/D");
  RecoTree->Branch("firstUnixTime",&firstUnixTime,"firstUnixTime/D");
  RecoTree->Branch("timeStamp",&timeStamp,"timeStamp/D");  
  RecoTree->Branch("runNum",&runNum,"runNum/I");
  RecoTree->Branch("Duration",&Duration,"Duration/D");
  RecoTree->Branch("FinalMinValue",&FinalMinValue,"FinalMinValue/D");
  RecoTree->Branch("Iterations",&Iterations,"Iterations/I");
  RecoTree->Branch("FinalTxCor_XYZ",FinalTxCor_XYZ,"FinalTxCor_XYZ[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR",FinalTxCor_ThPhR,"FinalTxCor_ThPhR[3]/D");
  RecoTree->Branch("InitialTxCor_XYZ",InitialTxCor_XYZ,"InitialTxCor_XYZ[3]/D");
  RecoTree->Branch("InitialTxCor_ThPhR",InitialTxCor_ThPhR,"InitialTxCor_ThPhR[3]/D");  
  RecoTree->Branch("dX",&dX,"dX/D");
  RecoTree->Branch("dY",&dY,"dY/D");
  RecoTree->Branch("dZ",&dZ,"dZ/D");
  RecoTree->Branch("dTh",&dTh,"dTh/D");
  RecoTree->Branch("dPh",&dPh,"dPh/D");
  RecoTree->Branch("dR",&dR,"dR/D");
  
  ceventTree->GetEntry(5);
  firstUnixTime=rawAtriEvPtr->unixTime;
  cout<<"firstUnixTime "<<setprecision(10)<<firstUnixTime<<endl;
  
  ceventTree->GetEntry(Event);
  cout<<Event<<" "<<rawAtriEvPtr->eventNumber<<endl;
  /////Right now only look at RF events which are not calpulser and software triggers
  if(rawAtriEvPtr->isCalpulserEvent()==false && rawAtriEvPtr->timeStamp>1000 && rawAtriEvPtr->isSoftwareTrigger()==false){
    ///Initialise the Useful event pointer to load the ARA event waveforms with the right calibration
    UsefulAtriStationEvent *realAtriEvPtr=new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
      
    eventNum=rawAtriEvPtr->eventNumber;
    unixTime=rawAtriEvPtr->unixTime;
    timeStamp=rawAtriEvPtr->timeStamp;   
    runNum=Run;
    
    //cout<<ievt<<" "<<V2nSum/(MCH-NumCutCh)<<" "<<CorScoreSum/(MCH-NumCutCh)<<endl;
    double ChHitTime[2][TotalAntennasRx];
    int IgnoreCh[2][TotalAntennasRx];
    double ChSNR[2][TotalAntennasRx];	
    double ChAmp[2][TotalAntennasRx];

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


    for(int ich=0;ich<TotalAntennasRx;ich++){
      Bool_t SkipCh=false;
      for(int icut=0; icut<NumCutCh; icut++){
	if(ich==CutCh[icut]){
	  SkipCh=true;
	}
      }
      for(int iray=0;iray<2;iray++){
	IgnoreCh[iray][ich]=1;
	if(ChHitTime[iray][ich]==0 || SkipCh==true){
	  IgnoreCh[iray][ich]=0;
	}
      }
    }

    UnixTimeSelected.push_back(unixTime-firstUnixTime);
    for(int ich=0; ich<MCH; ich++){
      if(IgnoreCh[0][ich]!=0 && IgnoreCh[1][ich]!=0){
	PwrSNR[0][ich].push_back(ChSNR[0][ich]);
	PwrSNR[1][ich].push_back(ChSNR[1][ich]);
	dtDR_Ch[ich].push_back(fabs(ChHitTime[0][ich]-ChHitTime[1][ich]));
	DRAmpRatio_Ch[ich].push_back(ChAmp[0][ich]/ChAmp[1][ich]);
      }else{
	PwrSNR[0][ich].push_back(ChSNR[0][ich]);  
	PwrSNR[1][ich].push_back(0);
	dtDR_Ch[ich].push_back(0);
	DRAmpRatio_Ch[ich].push_back(0);
      }
    }    

    double SPICE_Depth = gsl_spline_eval(spline_steffen, unixTime, spline_acc);
    InitialTxCor_XYZ[0]=TrueX;
    InitialTxCor_XYZ[1]=TrueY;
    InitialTxCor_XYZ[2]=-SPICE_Depth;

    FinalTxCor_XYZ[0]=0;
    FinalTxCor_XYZ[1]=0;
    FinalTxCor_XYZ[2]=0;
    
    DoInterferometery(InitialTxCor_XYZ, FinalTxCor_XYZ, InitialTxCor_XYZ, ExpectedPositionUncertainty, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, Duration, Iterations);      
      
    dX=InitialTxCor_XYZ[0]-FinalTxCor_XYZ[0];
    dY=InitialTxCor_XYZ[1]-FinalTxCor_XYZ[1];
    dZ=InitialTxCor_XYZ[2]-FinalTxCor_XYZ[2];
    
    InitialTxCor_ThPhR[0]=0;
    InitialTxCor_ThPhR[1]=0;
    InitialTxCor_ThPhR[2]=0;
    Interferometer::XYZtoThPhR(InitialTxCor_XYZ,InitialTxCor_ThPhR);
    InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
    InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 

    FinalTxCor_ThPhR[0]=0;
    FinalTxCor_ThPhR[1]=0;
    FinalTxCor_ThPhR[2]=0;
    Interferometer::XYZtoThPhR(FinalTxCor_XYZ, FinalTxCor_ThPhR);
    FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(180./Interferometer::pi);
    FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(180./Interferometer::pi); 

    dTh=InitialTxCor_ThPhR[0]-FinalTxCor_ThPhR[0];
    dPh=InitialTxCor_ThPhR[1]-FinalTxCor_ThPhR[1];
    dR=InitialTxCor_ThPhR[2]-FinalTxCor_ThPhR[2];
      
    cout<<"Final Reco Results are: dX="<<dX<<" ,dY="<<dY<<" ,dZ="<<dZ<<" |  Xtrue="<<InitialTxCor_XYZ[0]<<" ,Ytrue="<<InitialTxCor_XYZ[1]<<" ,Ztrue="<<InitialTxCor_XYZ[2]<<" | Xreco="<<FinalTxCor_XYZ[0]<<" ,Yreco="<<FinalTxCor_XYZ[1]<<" ,Zreco="<<FinalTxCor_XYZ[2]<<endl;
    cout<<"Final Reco Results are: dTh="<<dTh<<" ,dPh="<<dPh<<" ,dR="<<dR<<" |  Thtrue="<<InitialTxCor_ThPhR[0]<<" ,Phtrue="<<InitialTxCor_ThPhR[1]<<" ,Rtrue="<<InitialTxCor_ThPhR[2]<<" | Threco="<<FinalTxCor_ThPhR[0]<<" ,Phreco="<<FinalTxCor_ThPhR[1]<<" ,Rreco="<<FinalTxCor_ThPhR[2]<<endl;
    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;

     RecoTree->Fill();
    
    ///Delete the event pointer so we get a new pointer for the next event and there is no memory
    delete realAtriEvPtr;
      
  }////if condition for the single event  

  delete InputFile;

  OutputFile->Write();
  OutputFile->Close();
  
}
