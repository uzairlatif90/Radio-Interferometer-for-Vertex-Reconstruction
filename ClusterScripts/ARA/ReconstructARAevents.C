const int MCH=16;

#include "FFTtools.h"
#include "../../RET_int/Interferometer/Interferometer.cc"
#include <gsl/gsl_math.h>

TGraph *gCPtemp[2][16];

void ReadCPTemp(){

  {
    TString filename="../MyCorrCode/";
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
    TString filename="../MyCorrCode/";
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

void GetCPCor(Int_t iStation, vector <double> CalPulCor[3]){

  ///////////////////////////////////////////////////////////////
  Double_t antLocD5[2][3];
  Double_t antLocD6[2][3];
  
  AraStationInfo *stationInfo=new AraStationInfo(iStation,2019);
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

//this function calculates the SNR, signal amplitudes and noise RMSs of the voltage waveforms
Double_t getmyWaveformSNR(TGraph *gr){

  Int_t nBins = gr->GetN();
  Double_t *yVals = gr->GetY();
  
  double RMS=GetNoiseRMS(gr);
  
  ///find the signal amplitude
  Int_t trending=3;
  Double_t p2p=0;
  Int_t firstBin=0;
  Double_t y;

  Double_t newp2p=0;
  Int_t maxbin=0;
  for(Int_t i=0;i<nBins;i++){
    y=yVals[i];
    //cout<<i<<" "<<yVals[i]<<endl;
    if(i>0){
      if(y<yVals[i-1] && trending==0){
        if(TMath::Abs(y-yVals[firstBin]>p2p)){
          p2p=TMath::Abs(y-yVals[firstBin]);
        }
      }
      else if(y<yVals[i-1] && (trending==1 || trending==2)){
        trending=0;
        firstBin=i-1;
        if(TMath::Abs(y-yVals[firstBin]>p2p)){
          p2p=TMath::Abs(y-yVals[firstBin]);
        }
      }
      else if(y>yVals[i-1] && (trending==0 || trending==2)){
        trending=1;
        firstBin=i-1;
        if(TMath::Abs(y-yVals[firstBin]>p2p)){
          p2p=TMath::Abs(y-yVals[firstBin]);
	}
      }
      else if(y>yVals[i-1] && trending==1){
        if(TMath::Abs(y-yVals[firstBin]>p2p)){
          p2p=TMath::Abs(y-yVals[firstBin]);
	}
      }
      else if(y==yVals[i-1]){
	trending=2;
      }
      else if(trending==3){
        if(y<yVals[i-1]){
          trending=0;
          firstBin=0;
        }
        if(y>yVals[i-1]){
          trending=1;
          firstBin=0;
	}
      }
      else{
	std::cout << "trending cock up!" << std::endl;
	std::cout << "y " << y << " yVals[i] " << yVals[i] << " yVals[i-1] " << yVals[i-1] << std::endl;
        //return -1;
      }
      // if(p2p>newp2p){
      //   //cout<<"The value of firstbin is "<<firstBin<<" "<<p2p<<endl;
      // newp2p=p2p;
      //   maxbin=firstBin;
      // }
    }
  }
  
  //cout<<"The value of firstbin is "<<firstBin<<endl;
  p2p/=2.;
  Double_t SNR = p2p/RMS;

  //delete []yVals;
  
  return SNR;
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
    cout<<"We have two peaks "<<endl;

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
    cout<<"We have one peak "<<endl;

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


void ReconstructARAevents(Int_t StationId, char const *InputFileName, int Run, int Event){
// //void ReconstructARAevents(){
//   Int_t StationId=2;
//   char const *InputFileName="event9129.root";
//   int Run=9129;
//   int Event=1301;

  TString OutputFileName="";
  OutputFileName+="Run";
  OutputFileName+=Run;
  OutputFileName+="Event";
  OutputFileName+=Event;
  OutputFileName+=".root";   
  
  double ExpectedPositionUncertainty=5;//in m
  
  DeclareAntennaConfigARA(StationId);
  ReadCPTemp();
  
  ////initialise the event pointer to take data from the event tree
  RawAtriStationEvent *rawAtriEvPtr=0;

  ////initialise the AraEventCalibrator class  
  AraEventCalibrator *calib = AraEventCalibrator::Instance();
  
  ///Open the ARA root file that contains the data
  TFile *InputFile=new TFile(InputFileName);

  ///Create the output root file
  TFile *OutputFile=new TFile(OutputFileName,"RECREATE");
  
  ///Open the Tree in the file that contains data from all of the events
  TTree *ceventTree = (TTree*)InputFile->Get("eventTree");

  int runNumber;
  ///Set the Branch addresses
  ceventTree->SetBranchAddress("event",&rawAtriEvPtr);
  ceventTree->SetBranchAddress("run",&runNumber);
  cout<<"Opened the file and set the Branches"<<endl; 
  
  ///Get the total number of entries within the tree
  int nument=ceventTree->GetEntries();

  vector <double> PwrSNR[2][MCH];

  TGraph *gr[MCH];
  TGraph *grPwrEnv[MCH];
  TGraph *grCor[MCH];  

  double SNRV[8];
  double SNRH[8];
  
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

  TTree *RecoTree = new TTree("RecoTree","Reco info about Event");
  RecoTree->Branch("eventNum",&eventNum,"eventNum/I");
  RecoTree->Branch("unixTime",&unixTime,"unixTime/D");
  RecoTree->Branch("firstUnixTime",&firstUnixTime,"firstUnixTime/D");
  RecoTree->Branch("timeStamp",&timeStamp,"timeStamp/D");  
  RecoTree->Branch("runNum",&runNum,"runNum/I");
  RecoTree->Branch("Duration",&Duration,"Duration/D");
  RecoTree->Branch("FinalMinValue",&FinalMinValue,"FinalMinValue/D");
  RecoTree->Branch("Iterations",&Iterations,"Iterations/I");
  RecoTree->Branch("isCalpulserTrig",&isCalpulserTrig,"isCalpulserTrig/B");
  RecoTree->Branch("isSoftwareTrig",&isSoftwareTrig,"isSoftwareTrig/B");
  RecoTree->Branch("FinalTxCor_XYZ",FinalTxCor_XYZ,"FinalTxCor_XYZ[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR",FinalTxCor_ThPhR,"FinalTxCor_ThPhR[3]/D");
  RecoTree->Branch("InitialTxCor_XYZ",InitialTxCor_XYZ,"InitialTxCor_XYZ[3]/D");
  RecoTree->Branch("InitialTxCor_ThPhR",InitialTxCor_ThPhR,"InitialTxCor_ThPhR[3]/D");  
  RecoTree->Branch("VoltageSNR",VoltageSNR,"VoltageSNR[16]/D"); 
  // RecoTree->Branch("VoltageSNRV",VoltageSNRV,"VoltageSNRV[3]/D");
  // RecoTree->Branch("VoltageSNRV_Ch",VoltageSNRV_Ch,"VoltageSNRV_Ch[3]/I");
  // RecoTree->Branch("VoltageSNRH",VoltageSNRH,"VoltageSNRH[3]/D");
  // RecoTree->Branch("VoltageSNRH_Ch",VoltageSNRH_Ch,"VoltageSNRH_Ch[3]/I");
  RecoTree->Branch("CorScore",CorScore,"CorScore[16]/D");
  RecoTree->Branch("PowerSNR",PowerSNR,"PowerSNR[16]/D");
  
  ceventTree->GetEntry(5);
  firstUnixTime=rawAtriEvPtr->unixTime;
  cout<<"firstUnixTime "<<setprecision(10)<<firstUnixTime<<endl;
  
  ceventTree->GetEntry(Event);
  cout<<Event<<" "<<rawAtriEvPtr->eventNumber<<endl;

  ///Initialise the Useful event pointer to load the ARA event waveforms with the right calibration
  UsefulAtriStationEvent *realAtriEvPtr=new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
  isCalpulserTrig=rawAtriEvPtr->isCalpulserEvent();
  isSoftwareTrig=rawAtriEvPtr->isSoftwareTrigger();
    
  eventNum=rawAtriEvPtr->eventNumber;
  unixTime=rawAtriEvPtr->unixTime;
  timeStamp=rawAtriEvPtr->timeStamp;   
  runNum=Run;
  
  double ChHitTime[2][TotalAntennasRx];
  int IgnoreCh[2][TotalAntennasRx];
  double ChSNR[2][TotalAntennasRx];	
  double ChAmp[2][TotalAntennasRx];

  double dtDR_Avg=0;
  int IgnorePeakCount=0;

  double CutCh[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int NumChAvailable=0;
  for(int ich=0; ich<MCH; ich++){
    ///Get the Waveform from the data file for each channel
    TGraph *grdum=realAtriEvPtr->getGraphFromRFChan(ich);

    ///Interpolate the waveforms to ensure equal spacing between the samples
    TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.6);
    delete grdum;

    VoltageSNR[ich]=getmyWaveformSNR(gr);
    if(VoltageSNR[ich]<6){
      //cout<<ich<<" channel cut "<<VoltageSNR[ich]<<endl;
      CutCh[ich]=1;
    }else{
      NumChAvailable++;
    }
    
    delete gr;
  }
  CutCh[15]=1;
  
  for(int ich=0; ich<MCH; ich++){
    CorScore[ich]=0;
    
    if(CutCh[ich]==-1){
      ///Get the Waveform from the data file for each channel
      TGraph *grdum=realAtriEvPtr->getGraphFromRFChan(ich);

      ///Interpolate the waveforms to ensure equal spacing between the samples
      TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.6);
      delete grdum;

      //VoltageSNR[ich]=getmyWaveformSNR(gr);
	
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
	cout<<ich<<" "<<Dtime<<" "<<VoltageSNR[ich]<<" "<<DSNR<<endl;
      }

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
      CorScore[ich]=MaxCorScore;
	  
      delete gr;
      delete grPeakPoint;
    }

  }////channel loop

  ///separate out the V and H SNRs for voltage 
  for(Int_t cho=0; cho<8; cho++){
    SNRV[cho]= VoltageSNR[cho];
  }     
  for(Int_t ch=8; ch<16; ch++){
    SNRH[ch-8]=VoltageSNR[ch];
  }
    
  int nr=0;
  while(nr<3){
    VoltageSNRV[nr]=TMath::MaxElement(8,SNRV);
    VoltageSNRV_Ch[nr]=TMath::LocMax(8,SNRV);
    const Int_t dummy1=VoltageSNRV_Ch[nr];
    SNRV[dummy1]=0;
	    
    VoltageSNRH[nr]=TMath::MaxElement(8,SNRH);
    VoltageSNRH_Ch[nr]=TMath::LocMax(8,SNRH);
    const Int_t dummy2=VoltageSNRH_Ch[nr];
    SNRH[dummy2]=0;
    nr++;
  }

  for(int ich=0;ich<TotalAntennasRx;ich++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][ich]=1;
      if(ChHitTime[iray][ich]==0 || CutCh[ich]==1){
	IgnoreCh[iray][ich]=0;
      }
    }
  }

  for(int ich=0; ich<MCH; ich++){
    if(IgnoreCh[0][ich]!=0 && IgnoreCh[1][ich]!=0){
      PwrSNR[0][ich].push_back(ChSNR[0][ich]);
      PwrSNR[1][ich].push_back(ChSNR[1][ich]);
    }else{
      PwrSNR[0][ich].push_back(ChSNR[0][ich]);  
      PwrSNR[1][ich].push_back(0);
    }
  }

  if(NumChAvailable>=3){

    double GuessResultCor[3][3]; 
   
    if(rawAtriEvPtr->isCalpulserEvent()==true){
      vector <double> CalPulCor[3];
      GetCPCor(StationId, CalPulCor);
      Interferometer::GetApproximateMinUserCor(CalPulCor,GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
    }else{
      Interferometer::GetApproximateMinThPhR(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
      Interferometer::GetApproximateDistance(GuessResultCor,ExpectedPositionUncertainty,ChHitTime,IgnoreCh,ChSNR);
    }
    
    InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
    InitialTxCor_ThPhR[2]=GuessResultCor[0][2];
 
    Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ExpectedPositionUncertainty, ChHitTime, IgnoreCh, ChSNR, FinalMinValue, Duration, Iterations); 

    InitialTxCor_XYZ[0]=0;
    InitialTxCor_XYZ[1]=0;
    InitialTxCor_XYZ[2]=0;
    Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
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
  
    cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<endl;
    cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
    cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<Duration<<endl;
  }else{

    FinalTxCor_XYZ[0]=0;
    FinalTxCor_XYZ[1]=0;
    FinalTxCor_XYZ[2]=0;

    InitialTxCor_XYZ[0]=0;
    InitialTxCor_XYZ[1]=0;
    InitialTxCor_XYZ[2]=0;

    FinalTxCor_ThPhR[0]=0;
    FinalTxCor_ThPhR[1]=0;
    FinalTxCor_ThPhR[2]=0;

    InitialTxCor_ThPhR[0]=0;
    InitialTxCor_ThPhR[1]=0;
    InitialTxCor_ThPhR[2]=0;

    cout<<" Number of hit channels is less than 3 so no reconstruction was performed!!!"<<endl;
  }

  RecoTree->Fill();
    
  ///Delete the event pointer so we get a new pointer for the next event and there is no memory
  delete realAtriEvPtr;
  delete InputFile;

  OutputFile->Write();
  OutputFile->Close(); 
  
}
