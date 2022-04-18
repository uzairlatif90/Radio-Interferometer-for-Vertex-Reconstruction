const int MCH=16;

#include "FFTtools.h"
#include "../../Interferometer/Interferometer.cc"
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

TGraph *gCPtemp[2][16];


void ReadCPTemp(){

  {
    TString filename="../Interferometer/";
    filename+="CP_D6VPol_A2.root";
    TFile *f = TFile::Open(filename, "READ");
    for(int ich=0;ich<MCH;ich++){
      TString GetGraph="gr";
      GetGraph+=ich;
      f->GetObject(GetGraph, gCPtemp[0][ich]);

      int N=gCPtemp[0][ich]->GetN();
      int DummyBin=TMath::LocMax(N,gCPtemp[0][ich]->GetY());
      double xp,yp,x,y;
      double maxVoltage=gCPtemp[0][ich]->GetPoint(DummyBin,xp,yp);
      for(int ip=0;ip<N;ip++){
        gCPtemp[0][ich]->GetPoint(ip,x,y);
        gCPtemp[0][ich]->SetPoint(ip,x-xp,y);
      }
      
    }
    delete f;
  }
  {
    TString filename="../Interferometer/";
    filename+="CP_D6HPol_A2.root";
    TFile *f = TFile::Open(filename, "READ");
    for(int ich=0;ich<MCH;ich++){
      TString GetGraph="gr";
      GetGraph+=ich;
      f->GetObject(GetGraph, gCPtemp[1][ich]);

      int N=gCPtemp[0][ich]->GetN();
      int DummyBin=TMath::LocMax(N,gCPtemp[0][ich]->GetY());
      double xp,yp,x,y;
      double maxVoltage=gCPtemp[0][ich]->GetPoint(DummyBin,xp,yp);
      for(int ip=0;ip<N;ip++){
        gCPtemp[0][ich]->GetPoint(ip,x,y);
        gCPtemp[0][ich]->SetPoint(ip,x-xp,y);
      }
      
    }
    delete f;
  }

}

double integrateBinPower( TGraph *plot, int numBinsToIntegrate, vector<double> &integratedBins)
{
  int nPoints = plot->GetN();
  if (nPoints < numBinsToIntegrate){
    return 0;
  }

  //  Double_t *xVals = plot->GetX();                                                                                                          
  Double_t *yVals = plot->GetY();
  std::deque<double> integrator;
  double sum = 0.;
  for (int i = 0; i < numBinsToIntegrate; i++){
    integrator.push_back(pow(yVals[i], 2));
    sum = accumulate(integrator.begin(), integrator.end(), 0);
  }
  double max = 0.;
  integratedBins.push_back(sum);

  for (int i = 0+numBinsToIntegrate; i < nPoints; i++){

    sum = sum - integrator[0];
    integrator.pop_front();
    integrator.push_back(pow(yVals[i], 2));
    sum = sum + integrator[numBinsToIntegrate-1];

    integratedBins.push_back(sum);

    if (sum > max){
      max = sum;
    }
  }

  return max;
}

TGraph* makeIntegratedBinPowerGraphs(TGraph* graphsInput, int numBinsToIntegrate, string titles){
  
  vector<double> integratedBins;
  double maxIntPower = integrateBinPower(graphsInput, numBinsToIntegrate, integratedBins);
  double *volts = &integratedBins[0];
  TGraph* graphsOutput = new TGraph(integratedBins.size(), graphsInput->GetX(), volts);
  graphsOutput->GetXaxis()->SetTitle("Time (ns)");
  graphsOutput->GetYaxis()->SetTitle("Integrated Power (arb units)");
  graphsOutput->SetTitle("ich");
  integratedBins.clear();
  //    delete gIntPower;

  return graphsOutput;
}

void getAbsMaximum_N(TGraph* plot, int nPeaks, double timeApart, vector<double> &xs, vector<double> &ys)
{
  // xs.clear();
  // ys.clear();

  int nPoints = plot->GetN();
  if (nPoints < nPeaks){
    cerr << "Number of points in waveform is fewer than the number of peaks requested. Decreasing number of peaks requested to match number points." << endl;
    nPeaks = nPoints;
  } 

  double x_temp, y_temp;
  double y_good, x_good;
  int test;
  double y_upper;

  for (int iPeak = 0; iPeak < nPeaks; iPeak++){
    y_good = -9.0E99;
    y_upper = 9.0E99;
    if (iPeak > 0){
      y_upper = ys[iPeak-1];
    }
    for (int i = 0; i < nPoints; i++){
      test = plot->GetPoint(i, x_temp, y_temp);
      if (iPeak == 0){
	if (y_temp > y_good){
	  x_good = x_temp;
	  y_good = y_temp;
	}
      } 
      else {
	for (int iTimeTest = 0; iTimeTest < xs.size(); iTimeTest++){
	  if (y_temp > y_good && y_temp < y_upper && abs(x_temp - xs[iTimeTest]) > timeApart){
	    x_good = x_temp;
	    y_good = y_temp;
	  }
	}
      }
    }
    xs.push_back(x_good);
    ys.push_back(y_good);
    //cout << iPeak << " : " << xs[iPeak] << " : " << ys[iPeak] << endl;
  }
  //return;
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

  // for(int ipeak=0;ipeak<2;ipeak++){
  //   vector <double> RefinePeak[3];
  //   for(int isample=PowerPeakBin[ipeak]-DeltaSamples; isample<=PowerPeakBin[ipeak]+DeltaSamples; isample++){
  //     double xp,yp;
  //     grPwrEnvOriginal->GetPoint(isample,xp,yp);
  //     RefinePeak[0].push_back(xp);
  //     RefinePeak[1].push_back(yp);
  //     RefinePeak[2].push_back(isample);
  //   }
      
  //   DummyBin=TMath::LocMax(RefinePeak[1].size(),RefinePeak[1].data());
  //   PowerPeakAmp[ipeak]=RefinePeak[1][DummyBin];
  //   PowerPeakTime[ipeak]=RefinePeak[0][DummyBin];
  //   PowerPeakBin[ipeak]=RefinePeak[2][DummyBin];
  //   RefinePeak[0].clear();
  //   RefinePeak[1].clear();
  //   RefinePeak[2].clear();
  // }
  

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
      
      // vector <double> RefinePeak[3];
      // for(int isample=PowerPeakBin[ipeak]-20; isample<=PowerPeakBin[ipeak]+20; isample++){
      // 	if(isample>-1){
      // 	  double xp,yp;
      // 	  grPwrEnvOriginal->GetPoint(isample,xp,yp);
      // 	  if(isample<PowerPeakBin[ipeak]-5 || isample>PowerPeakBin[ipeak]+5){
      // 	    RefinePeak[0].push_back(xp);
      // 	    RefinePeak[1].push_back(yp);
      // 	    RefinePeak[2].push_back(isample);
      // 	  }
      // 	}
      // }
      
      // DummyBin=TMath::LocMax(RefinePeak[1].size(),RefinePeak[1].data());
      // PowerPeakAmp[ipeak+2]=RefinePeak[1][DummyBin];
      // PowerPeakTime[ipeak+2]=RefinePeak[0][DummyBin];
      // PowerPeakBin[ipeak+2]=RefinePeak[2][DummyBin];
      // RefinePeak[0].clear();
      // RefinePeak[1].clear();
      // RefinePeak[2].clear(); 
    
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
     
      // vector <double> RefinePeak[3];
      // for(int isample=PowerPeakBin[ipeak]-20; isample<=PowerPeakBin[ipeak]+20; isample++){
      // 	if(isample>-1){
      // 	  double xp,yp;
      // 	  grPwrEnvOriginal->GetPoint(isample,xp,yp);
      // 	  if(isample<PowerPeakBin[ipeak]-5 || isample>PowerPeakBin[ipeak]+5){
      // 	    RefinePeak[0].push_back(xp);
      // 	    RefinePeak[1].push_back(yp);
      // 	    RefinePeak[2].push_back(isample);
      // 	  }
      // 	}
      // }
      
      // DummyBin=TMath::LocMax(RefinePeak[1].size(),RefinePeak[1].data());
      // PowerPeakAmp[ipeak+2]=RefinePeak[1][DummyBin];
      // PowerPeakTime[ipeak+2]=RefinePeak[0][DummyBin];
      // PowerPeakBin[ipeak+2]=RefinePeak[2][DummyBin];
      // RefinePeak[0].clear();
      // RefinePeak[1].clear();
      // RefinePeak[2].clear(); 
    
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

void ReconstructSPICEevents(int StationId,char const *InputFileName, int Run, int Event){

  DeclareAntennaConfigARA(StationId);
  LoadDepthFile();
  ReadCPTemp();
  
  TString OutputFileName="./output/";
  OutputFileName+="Run";
  OutputFileName+=Run;
  OutputFileName+="Event";
  OutputFileName+=Event;
  OutputFileName+=".root";

  //A2
  double TrueX=-456.721;
  double TrueY=-2353;
  //A5
  // double TrueX=348.963;
  // double TrueY=-4152.76;

  double ExpectedPositionUncertainty=5;//in m
  double GetActualInitialCondition=true;
  bool RefineRecoResults=false;
  
  ////initialise the event pointer to take data from the event tree
  RawAtriStationEvent *rawAtriEvPtr=0;

  ////initialise the AraEventCalibrator class  
  //AraEventCalibrator *calib = AraEventCalibrator::Instance();
  
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

  vector <double> UnixTime[MCH];
  vector <double> V2n[MCH];
  vector <double> PwrSNR[2][MCH];
  
  vector <double> DRAmpRatio_Ch[MCH];
  vector <double> SPICEReco[3];
  vector <double> UnixTimeSelected;
  vector <double> UnixTimeReco;
  vector <double> DurationReco;

  vector <double> V2nAvg;
  vector <double> CorScoreAvg;

  TGraph *grWF[MCH];
  
  TGraph *gr[MCH];
  TGraph *grPwrEnv[MCH];
  TGraph *grCor[MCH];  

  double SNRV[8];
  double SNRH[8];
  
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
  int IsItBelowStation;
  double ArrivalDirection[3];
 
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

  TTree *RecoTree = new TTree("RecoTree","Reco info about Event");
  RecoTree->Branch("eventNum",&eventNum,"eventNum/D");
  RecoTree->Branch("unixTime",&unixTime,"unixTime/D");
  RecoTree->Branch("firstUnixTime",&firstUnixTime,"firstUnixTime/D");
  RecoTree->Branch("timeStamp",&timeStamp,"timeStamp/D");  
  RecoTree->Branch("runNum",&runNum,"runNum/I");
  RecoTree->Branch("DurationTotal",&DurationTotal,"DurationTotal/D");
  RecoTree->Branch("DurationReconstruction",&DurationReconstruction,"DurationReconstruction/D");
  RecoTree->Branch("DurationInitialCondition",&DurationInitialCondition,"DurationInitialCondition/D");
  RecoTree->Branch("FinalMinValue",&FinalMinValue,"FinalMinValue/D");
  RecoTree->Branch("Iterations",&Iterations,"Iterations/I");
  RecoTree->Branch("isCalpulserTrig",&isCalpulserTrig,"isCalpulserTrig/O");
  RecoTree->Branch("isSoftwareTrig",&isSoftwareTrig,"isSoftwareTrig/O");

  RecoTree->Branch("FinalTxCor_XYZ",FinalTxCor_XYZ,"FinalTxCor_XYZ[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR",FinalTxCor_ThPhR,"FinalTxCor_ThPhR[3]/D");
  RecoTree->Branch("FinalTxCor_XYZ_fR",FinalTxCor_XYZ_fR,"FinalTxCor_XYZ_fR[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR_fR",FinalTxCor_ThPhR_fR,"FinalTxCor_ThPhR_fR[3]/D");
  RecoTree->Branch("InitialTxCor_XYZ",InitialTxCor_XYZ,"InitialTxCor_XYZ[3]/D");
  RecoTree->Branch("InitialTxCor_ThPhR",InitialTxCor_ThPhR,"InitialTxCor_ThPhR[3]/D");  
  RecoTree->Branch("VoltageSNR",VoltageSNR,"VoltageSNR[16]/D"); 
  // RecoTree->Branch("VoltageSNRV",VoltageSNRV,"VoltageSNRV[3]/D");
  // RecoTree->Branch("VoltageSNRV_Ch",VoltageSNRV_Ch,"VoltageSNRV_Ch[3]/I");
  // RecoTree->Branch("VoltageSNRH",VoltageSNRH,"VoltageSNRH[3]/D");
  // RecoTree->Branch("VoltageSNRH_Ch",VoltageSNRH_Ch,"VoltageSNRH_Ch[3]/I");
  RecoTree->Branch("CorScore",CorScore,"CorScore[16]/D");
  RecoTree->Branch("ChHitTime",ChHitTime,"ChHitTime[2][16]/D");
  RecoTree->Branch("IgnoreCh",IgnoreCh,"IgnoreCh[2][16]/I");  
  RecoTree->Branch("dXYZ",dXYZ,"dXYZ[3]/D");
  RecoTree->Branch("dThPhR",dThPhR,"dThPhR[3]/D");
  
  ceventTree->GetEntry(5);
  firstUnixTime=rawAtriEvPtr->unixTime;
  cout<<"firstUnixTime "<<setprecision(10)<<firstUnixTime<<endl;
  
  ceventTree->GetEntry(Event);
  cout<<Event<<" "<<rawAtriEvPtr->eventNumber<<endl;
  /////Right now only look at RF events which are not calpulser and software triggers
  if(rawAtriEvPtr->isCalpulserEvent()==false && rawAtriEvPtr->timeStamp>1000 && rawAtriEvPtr->isSoftwareTrigger()==false){
    ///Initialise the Useful event pointer to load the ARA event waveforms with the right calibration
    UsefulAtriStationEvent *realAtriEvPtr=new UsefulAtriStationEvent(rawAtriEvPtr, AraCalType::kLatestCalib);
    isCalpulserTrig=rawAtriEvPtr->isCalpulserEvent();
    isSoftwareTrig=rawAtriEvPtr->isSoftwareTrigger();

    eventNum=rawAtriEvPtr->eventNumber;
    unixTime=rawAtriEvPtr->unixTime;
    timeStamp=rawAtriEvPtr->timeStamp;   
    runNum=Run;
    
    double ChSNR[2][TotalAntennasRx];	
    double ChAmp[2][TotalAntennasRx];

    double CutCh[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double FaultyCh[16]={-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,1,-1,1};//A2
    //double FaultyCh[16]={-1,-1,-1,1,-1,-1,-1,-1,1,1,1,1,1,1,1,1};//A5
    int NumChAvailable=0;
    int NumChAvailableV=0;
    int NumChAvailableH=0;
     
    for(int ich=0; ich<MCH; ich++){

              
      //Get the Waveform from the data file for each channel
      TGraph *grdum=usefulAtriEvPtr->getGraphFromRFChan(ich);

      //Interpolate the waveforms to ensure equal spacing between the samples
      TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.01);
      delete grdum;

      //Initiate values to set up the peak finding algorithm
      int numBinsToIntegrate=10;
      string title="ich";
      TGraph* grIntPower = makeIntegratedBinPowerGraphs(gr, numBinsToIntegrate,  title);

      vector<double> vvHitTimes; // vector of vector of hit times (first index is antenna, second index is hit time)
      vector<double> vvPeakIntPowers; // vector of vector of peak values
     
      //Find the peaks
      int numSearchPeaks=2; //only look for two peaks
      double peakSeparation=30.0;
      getAbsMaximum_N(grIntPower, numSearchPeaks, peakSeparation, vvHitTimes, vvPeakIntPowers);
      
      //Calculate the noise RMS for the waveform
      double noiserms=GetNoiseRMS(grIntPower);

    
      //Fill in the peak hit times and peak SNR values
      ChHitTime[0][ich]=vvHitTimes[0];
      ChHitTime[1][ich]=vvHitTimes[1];
      ChSNR[0][ich]=vvPeakIntPowers[0]/noiserms;
      ChSNR[1][ich]=vvPeakIntPowers[1]/noiserms;
      if(ChSNR[1][ich]/ChSNR[0][ich]<0.2){
	ChSNR[1][ich]=10;
      }     
      //Ignore the second peak for now by giving it a below threshold SNR
      ChSNR[1][ich]=10;
    
      //Filter out channels which have low SNR or are faulty
      if( (ChSNR[0][ich]<25 && ChSNR[1][ich]<25) || FaultyCh[ich]==1){
  	//if(FaultyCh[ich]==1){
	//cout<<ich<<" channel cut "<<VoltageSNR[ich]<<endl;
	CutCh[ich]=1;
	VoltageSNR[ich]=0;
      }else{
	NumChAvailable++;
	if(ich<8){
	  NumChAvailableV++;
	}else{
	  NumChAvailableH++;
	}
      }
       
      delete gr;
    }
    
    for(int ich=0; ich<MCH; ich++){
      CorScore[ich]=0;
      
      if(CutCh[ich]==-1){
	// ///Get the Waveform from the data file for each channel
	// TGraph *grdum=realAtriEvPtr->getGraphFromRFChan(ich);

	// ///Interpolate the waveforms to ensure equal spacing between the samples
	// TGraph *gr=FFTtools::getInterpolatedGraph(grdum,0.6);
	// delete grdum;

	// grPwrEnv[ich]=FFTtools::getSimplePowerEnvelopeGraph(gr);
	// TGraph *grPeakPoint=new TGraph();
	
	// PeakFinder(grPwrEnv[ich], grPeakPoint);
	  
	// int NumOfPulses=grPeakPoint->GetN();
	// double Dtime=0,Rtime=0;
	// double DAmp=0,RAmp=0;
	// double DSNR=0,RSNR=0;
	// int IgnorePeakD=0,IgnorePeakR=0;
	// if(NumOfPulses==6){
	double Dtime=0,Rtime=0;
	double DSNR=0,RSNR=0;
	
	//If both peaks are present
	if(ChSNR[0][ich]>=25 && ChSNR[1][ich]>=25){

	  // double xpeak,ypeak;
	  // grPeakPoint->GetPoint(0,xpeak,ypeak);
	  // Dtime=xpeak;
	  // DAmp=ypeak;
	  // grPeakPoint->GetPoint(2,xpeak,ypeak);
	  // DSNR=ypeak;
	  // grPeakPoint->GetPoint(4,xpeak,ypeak);
	  // IgnorePeakD=xpeak;
	   
	  // grPeakPoint->GetPoint(1,xpeak,ypeak);
	  // Rtime=xpeak;
	  // RAmp=ypeak;
	  // grPeakPoint->GetPoint(3,xpeak,ypeak);
	  // RSNR=ypeak;
	  // grPeakPoint->GetPoint(5,xpeak,ypeak);
	  // IgnorePeakR=xpeak;
	    
	  // if(Dtime>Rtime){
	  //   swap(Dtime,Rtime);
	  //   swap(DAmp,RAmp);
	  //   swap(DSNR,RSNR);
	  //   swap(IgnorePeakD,IgnorePeakR);
	  // }

	  // if(IgnorePeakD==0){
	  //   Dtime=0;
	  //   DSNR=0;
	  //   DAmp=0;
	  //   IgnorePeakCount++;
	  // }

	  // if(IgnorePeakR==0){
	  //   Rtime=0;
	  //   RSNR=0;
	  //   RAmp=0;
	  //   IgnorePeakCount++;
	  // }
	      
	  // ChHitTime[0][ich]=Dtime;
	  // ChHitTime[1][ich]=Rtime;
	  // ChSNR[0][ich]=DSNR;
	  // ChSNR[1][ich]=RSNR;
	  // ChAmp[0][ich]=DAmp;
	  // ChAmp[1][ich]=RAmp;

	  Dtime=ChHitTime[0][ich];
	  Rtime=ChHitTime[1][ich];
	  DSNR=ChSNR[0][ich];
	  RSNR=ChSNR[1][ich];
	  
	  if(Dtime>Rtime){
	    swap(Dtime,Rtime);
	    swap(ChSNR[0][ich],ChSNR[1][ich]);
	  }
  
	  ChHitTime[0][ich]=Dtime;
	  ChHitTime[1][ich]=Rtime;

	  cout<<ich<<" two peak "<<Dtime<<" "<<Rtime<<" "<<ChSNR[0][ich]<<" "<<ChSNR[1][ich]<<endl;
	  //cout<<ich<<" two peaks "<<VoltageSNR[ich]<<endl;
	}else{
	  // double xpeak,ypeak;
	  // grPeakPoint->GetPoint(0,xpeak,ypeak);
	  // Dtime=xpeak;
	  // DAmp=ypeak;
	  // grPeakPoint->GetPoint(1,xpeak,ypeak);
	  // DSNR=ypeak;
	  // grPeakPoint->GetPoint(2,xpeak,ypeak);
	  // IgnorePeakD=xpeak;
	    
	  // if(IgnorePeakD==0){
	  //   Dtime=0;
	  //   DSNR=0;
	  //   DAmp=0;
	  //   IgnorePeakCount++;
	  // }
	  // ChHitTime[0][ich]=Dtime;
	  // ChHitTime[1][ich]=0;
	  // ChSNR[0][ich]=DSNR;
	  // ChSNR[1][ich]=0;
	  // ChAmp[0][ich]=DAmp;
	  // ChAmp[1][ich]=0;
	  // cout<<ich<<" one peak "<<VoltageSNR[ich]<<endl;
	  
	  Dtime=ChHitTime[0][ich];
	  
	  ChHitTime[0][ich]=Dtime;
	  ChHitTime[1][ich]=0;
	  ChSNR[1][ich]=0;
	  cout<<ich<<" one peak "<<Dtime<<" "<<ChSNR[0][ich]<<endl;
	}
	  
	int N=gr->GetN();
	int DummyBin=TMath::LocMax(N,gr->GetY());
	double xp,yp,x,y;
	double maxVoltage=gr->GetPoint(DummyBin,xp,yp);
	for(int ip=0;ip<N;ip++){
	  gr->GetPoint(ip,x,y);
	  gr->SetPoint(ip,x-xp,y);
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
      if(IgnoreCh[0][cho]==1){
	SNRV[cho]= VoltageSNR[cho];
      }else{
	SNRV[cho]= 0;
      }
    }     
    for(Int_t ch=8; ch<16; ch++){
      if(IgnoreCh[0][ch]==1){
	SNRH[ch-8]=VoltageSNR[ch];
      }else{
	SNRH[ch-8]=0;
      }
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
  
    Double_t NumSinglePeak=0;
    Double_t NumTotalChannels=0;
    for(int ich=0;ich<TotalAntennasRx;ich++){
      for(int iray=0;iray<2;iray++){
	IgnoreCh[iray][ich]=1;
	if(ChHitTime[iray][ich]==0 || CutCh[ich]==1){
	  IgnoreCh[iray][ich]=0;
	}
      
      }
      if(IgnoreCh[0][ich]==1){
	NumTotalChannels++;
      }
      if((IgnoreCh[0][ich]==1 && IgnoreCh[1][ich]==0)){
	NumSinglePeak++;
      }
    } 

    Double_t FractionSinglePeak=NumSinglePeak/NumTotalChannels;
    if(FractionSinglePeak>=0.85){
      for(int ich=0;ich<TotalAntennasRx;ich++){
        if(IgnoreCh[1][ich]==1){
	IgnoreCh[1][ich]=0;
	if(ChAmp[1][ich]>ChAmp[0][ich]){
	  swap(ChHitTime[0][ich], ChHitTime[1][ich]);
	  swap(ChSNR[0][ich], ChSNR[1][ich]);
	}
        }
      
      }
    }

    vector <double> ChHitTimev[2]; ////Channel Hit Time
    vector <int> IgnoreChv[2];
    vector <double> ChSNRv[2]; 
    ChHitTimev[0].resize(TotalAntennasRx);
    ChHitTimev[1].resize(TotalAntennasRx);
    IgnoreChv[0].resize(TotalAntennasRx);
    IgnoreChv[1].resize(TotalAntennasRx);
    ChSNRv[0].resize(TotalAntennasRx);
    ChSNRv[1].resize(TotalAntennasRx);

    for(int ich=0;ich<TotalAntennasRx;ich++){
      for(int iray=0;iray<2;iray++){
	ChHitTimev[iray][iRx]=ChHitTime[iray][iRx];
	IgnoreChv[iray][iRx]=IgnoreCh[iray][iRx];
	ChSNRv[iray][iRx]=ChSNR[iray][iRx];
      }
    }
    
    if(NumChAvailable>=4){
      double GuessResultCor[3][3]; 
      double MinimizerRadialWidth;
        
      IsItBelowStation=Interferometer::IsItAboveOrBelow(ChHitTimev,IgnoreChv) ;
      
      Interferometer::GetRecieveAngle(ChHitTimev, IgnoreChv, ChSNRv, IsItBelowStation, 100, ArrivalDirection);
      
      bool CheckTrigger=true;
 
      auto t1 = std::chrono::high_resolution_clock::now();
      double SPICE_Depth = gsl_spline_eval(spline_steffen, unixTime, spline_acc);
 
      if(GetActualInitialCondition==true){  

	vector <double> timeRay[2];
	vector <int> IgnoreChB[2];
	
	InitialTxCor_XYZ[0]=TrueX -AvgAntennaCoordRx[0];
	InitialTxCor_XYZ[1]=TrueY -AvgAntennaCoordRx[1];
	InitialTxCor_XYZ[2]=-SPICE_Depth -AvgAntennaCoordRx[2];;
    
	Interferometer::XYZtoThPhR(InitialTxCor_XYZ,InitialTxCor_ThPhR);
	InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
	InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 

	GuessResultCor[0][0]=InitialTxCor_ThPhR[0];
	GuessResultCor[0][1]=InitialTxCor_ThPhR[1];
	GuessResultCor[0][2]=InitialTxCor_ThPhR[2];
      
	MinimizerRadialWidth=100;

	Interferometer::GenerateChHitTimeAndCheckHits(InitialTxCor_XYZ,timeRay,IgnoreChB);
	CheckTrigger=Interferometer::CheckTrigger(IgnoreChB);

      }else{
	Interferometer::GetApproximateMinThPhR(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,10);      
	Interferometer::GetApproximateDistance(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,10);
	MinimizerRadialWidth=500;
      }
      auto t2 = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
      DurationInitialCondition=duration/1000;

      InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
      InitialTxCor_ThPhR[2]=GuessResultCor[0][2];
      
      if(CheckTrigger==true){
 
	Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ExpectedPositionUncertainty, ChHitTimev, IgnoreChv, ChSNRv, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth); 
	
	if(RefineRecoResults==true){
	  cout<<"1st try Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]*(180./Interferometer::pi)<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]*(180./Interferometer::pi)<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
	  
	  InitialTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
	  InitialTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180);
	  InitialTxCor_ThPhR[2]=FinalTxCor_ThPhR[2];  
	  
	  Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTimev, IgnoreChv, ChSNRv, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50); 
	}
	
	double FixedR=sqrt(SPICE_Depth*SPICE_Depth+ TrueX*TrueX + TrueY*TrueY);
	Interferometer::GetRecoFixedR(GuessResultCor[0], FinalTxCor_ThPhR_fR,ChHitTimev, IgnoreChv, ChSNRv, FixedR);

	DurationTotal=DurationInitialCondition+DurationReconstruction;

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

	FinalTxCor_XYZ_fR[0]=0;
	FinalTxCor_XYZ_fR[1]=0;
	FinalTxCor_XYZ_fR[2]=0;
	FinalTxCor_ThPhR_fR[0]=FinalTxCor_ThPhR_fR[0]*(Interferometer::pi/180);
	FinalTxCor_ThPhR_fR[1]=FinalTxCor_ThPhR_fR[1]*(Interferometer::pi/180); 
	Interferometer::ThPhRtoXYZ(FinalTxCor_ThPhR_fR, FinalTxCor_XYZ_fR);
	FinalTxCor_ThPhR_fR[0]=FinalTxCor_ThPhR_fR[0]*(180./Interferometer::pi);
	FinalTxCor_ThPhR_fR[1]=FinalTxCor_ThPhR_fR[1]*(180./Interferometer::pi);

	for(int i=0;i<3;i++){
	  dXYZ[i]=InitialTxCor_XYZ[i]-FinalTxCor_XYZ[i]; 
	  dThPhR[i]=InitialTxCor_ThPhR[i]-FinalTxCor_ThPhR[i]; 
	}
    
	cout<<"Final Reco Results are: dX="<<dXYZ[0]<<" ,dY="<<dXYZ[1]<<" ,dZ="<<dXYZ[2]<<" |  Xtrue="<<InitialTxCor_XYZ[0]<<" ,Ytrue="<<InitialTxCor_XYZ[1]<<" ,Ztrue="<<InitialTxCor_XYZ[2]<<" | Xreco="<<FinalTxCor_XYZ[0]<<" ,Yreco="<<FinalTxCor_XYZ[1]<<" ,Zreco="<<FinalTxCor_XYZ[2]<<endl;
	cout<<"Final Reco Results are: dTh="<<dThPhR[0]<<" ,dPh="<<dThPhR[1]<<" ,dR="<<dThPhR[2]<<" |  Thtrue="<<InitialTxCor_ThPhR[0]<<" ,Phtrue="<<InitialTxCor_ThPhR[1]<<" ,Rtrue="<<InitialTxCor_ThPhR[2]<<" | Threco="<<FinalTxCor_ThPhR[0]<<" ,Phreco="<<FinalTxCor_ThPhR[1]<<" ,Rreco="<<FinalTxCor_ThPhR[2]<<endl;
	cout<<"Final Reco Results for fixed R are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<FixedR<<" | Th_reco="<<FinalTxCor_ThPhR_fR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR_fR[1]<<" ,R_reco="<<FixedR<<endl;
	cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<DurationTotal<<endl;
      }else{
	
	//If no reconstruction was performed then fill the final arrays with zeroes
	ArrivalDirection[0]=0;
	ArrivalDirection[1]=0;
	ArrivalDirection[2]=0;
	
	FinalTxCor_XYZ[0]=0;
	FinalTxCor_XYZ[1]=0;
	FinalTxCor_XYZ[2]=0;
	
	InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
	InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 
	
	FinalTxCor_ThPhR[0]=0;
	FinalTxCor_ThPhR[1]=0;
	FinalTxCor_ThPhR[2]=0;

	cout<<" Initial condition is in shadow zone no reconstruction was performed!!!"<<endl;
      }

    }else{

      
      //If no reconstruction was performed then fill the final arrays with zeroes
      ArrivalDirection[0]=0;
      ArrivalDirection[1]=0;
      ArrivalDirection[2]=0;
      
      FinalTxCor_XYZ[0]=0;
      FinalTxCor_XYZ[1]=0;
      FinalTxCor_XYZ[2]=0;

      InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
      InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 
      
      FinalTxCor_ThPhR[0]=0;
      FinalTxCor_ThPhR[1]=0;
      FinalTxCor_ThPhR[2]=0;

      cout<<" Number of hit channels is less than 3 so no reconstruction was performed!!!"<<endl;
    }  

    RecoTree->Fill();
    RecoTree->Write();
    
    ///Delete the event pointer so we get a new pointer for the next event and there is no memory
    delete realAtriEvPtr;
      
  }////if condition for the single event  

  delete InputFile;

  OutputFile->Write();
  OutputFile->Close();
  
}
