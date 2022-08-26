const int MCH=16;
UsefulAtriStationEvent *usefulAtriEvPtr;

// AraSim includes 
#include "Detector.h"
#include "Event.h"
#include "Position.h"

#include "FFTtools.h"
#include "Interferometer.cc"

//TGraph *grTemp[MCH];

int getRFChannel(int string, int channel){ //mapping from AraSim to RF channel chain
  int RFChannel = 0;
  if(string == 0){
    if(channel == 0){
      RFChannel = 5;
    }
    if(channel == 1){
      RFChannel = 13;
    }
    if(channel == 2){
      RFChannel = 1;
    }
    if(channel == 3){
      RFChannel = 9;
    }
  }   
  if(string == 1){
    if(channel == 0){
      RFChannel = 6;
    }
    if(channel == 1){
      RFChannel = 14;
    }
    if(channel == 2){
      RFChannel = 2;
    }
    if(channel == 3){
      RFChannel = 10;
    }
  }     
  if(string == 2){
    if(channel == 0){
      RFChannel = 7;
    }
    if(channel == 1){
      RFChannel = 15;
    }
    if(channel == 2){
      RFChannel = 3;
    }
    if(channel == 3){
      RFChannel = 11;
    }
  }     
  if(string == 3){
    if(channel == 0){
      RFChannel = 4;
    }
    if(channel == 1){
      RFChannel = 12;
    }
    if(channel == 2){
      RFChannel = 0;
    }
    if(channel == 3){
      RFChannel = 8;
    }
  }
    
  return RFChannel;
}

// void ReadTemplate(){

//   TString filename="/data/user/ulatif/AraSimReco/";
//   filename+="Template.root";
//   TFile *f = TFile::Open(filename, "READ");
//   for(int ich=0;ich<MCH;ich++){
//     TString GetGraph="GraphCh";
//     GetGraph+=ich;
//     f->GetObject(GetGraph, grTemp[ich]);
//   }
//   delete f;

// }


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

double *get_detector_cog(Detector *detectorPtr){
  double *average_position=new double[3];
  average_position[0]=0;
  average_position[1]=0;
  average_position[2]=0;

  int n_ant = 0;
  for(int s=0; s<4; s++){
    for(int a=0; a<4; a++){
     
      average_position[0]+=detectorPtr->stations[0].strings[s].antennas[a].GetX();
      average_position[1]+=detectorPtr->stations[0].strings[s].antennas[a].GetY();
      average_position[2]+=detectorPtr->stations[0].strings[s].antennas[a].GetZ();
      n_ant++;
    }
  }
  
  average_position[0]/=16;
  average_position[1]/=16;
  average_position[2]/=16;

  return average_position;
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

void ReconstructAraSimevents(int Run, char const *InputFileName,char const *OutputFileName){

  ///Open the ARA root file that contains the data
  TFile *InputFile = TFile::Open(InputFileName);

  // data like      
  int runNumber;
  double weight; 
 
  TTree *eventTree = (TTree*) InputFile->Get("eventTree");
  eventTree->SetBranchAddress("UsefulAtriStationEvent", &usefulAtriEvPtr);
  eventTree->SetBranchAddress("weight", &weight);

  ///Get the total number of entries within the tree
  int nument=eventTree->GetEntries();
  
  // simulation                                                 
  TTree *simSettingsTree;
  Detector *detector = 0;
  simSettingsTree=(TTree*) InputFile->Get("AraTree");
  simSettingsTree->SetBranchAddress("detector", &detector);
  simSettingsTree->GetEntry(0);

  TTree *simTree;
  Event *eventPtr = 0; // it is apparently incredibly important that this be initialized to zero.
  Report *reportPtr = 0; // it is apparently incredibly important that this be initialized to zero.
  simTree=(TTree*) InputFile->Get("AraTree2");
  simTree->SetBranchAddress("event", &eventPtr);
  simTree->SetBranchAddress("report", &reportPtr);

  ////Declare Antenna Geometry which is required by the interferometer
  DeclareAntennaConfigAraSim();
  //ReadTemplate();

  //Set the name for the output file
  // TString OutputFileName="/data/user/ulatif/AraSimReco/output/";
  // OutputFileName+="Run";
  // OutputFileName+=Run;
  // OutputFileName+=".root";

  ///Create the output root file
  TFile *OutputFile=new TFile(OutputFileName,"RECREATE");
  
  //create the varaibles that will be stored in the tree that will be stored in the output root file
  int eventNum;
  int runNum;
  
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
  double VoltageSNR[16]; 
  
  double ChSNR[2][16];	
  double ChHitTime[2][16];
  double ChHitTime_WF[2][16];
  int IgnoreCh[2][16];

  double neutrinoWeight;
  int IsItBelowStation;

  double pnu;

  double ArrivalDirection[3];
  
  double timeRay_AS[2][16];
  
  TTree *RecoTree = new TTree("RecoTree","Reco info about Event");
  RecoTree->Branch("eventNum",&eventNum,"eventNum/I");
  RecoTree->Branch("runNum",&runNum,"runNum/I");
  RecoTree->Branch("DurationTotal",&DurationTotal,"DurationTotal/D");
  RecoTree->Branch("DurationReconstruction",&DurationReconstruction,"DurationReconstruction/D");
  RecoTree->Branch("DurationInitialCondition",&DurationInitialCondition,"DurationInitialCondition/D");
  RecoTree->Branch("FinalMinValue",&FinalMinValue,"FinalMinValue/D");
  RecoTree->Branch("Iterations",&Iterations,"Iterations/I");
  RecoTree->Branch("FinalTxCor_XYZ",FinalTxCor_XYZ,"FinalTxCor_XYZ[3]/D");
  RecoTree->Branch("FinalTxCor_ThPhR",FinalTxCor_ThPhR,"FinalTxCor_ThPhR[3]/D");
  RecoTree->Branch("TrueTxCor_XYZ",TrueTxCor_XYZ,"TrueTxCor_XYZ[3]/D");
  RecoTree->Branch("TrueTxCor_ThPhR",TrueTxCor_ThPhR,"TrueTxCor_ThPhR[3]/D");
  RecoTree->Branch("InitialTxCor_XYZ",InitialTxCor_XYZ,"InitialTxCor_XYZ[3]/D");
  RecoTree->Branch("InitialTxCor_ThPhR",InitialTxCor_ThPhR,"InitialTxCor_ThPhR[3]/D");  
  RecoTree->Branch("ChSNR",ChSNR,"ChSNR[2][16]/D");
  RecoTree->Branch("ChHitTime",ChHitTime,"ChHitTime[2][16]/D");
  RecoTree->Branch("IgnoreCh",IgnoreCh,"IgnoreCh[2][16]/I"); 
  RecoTree->Branch("neutrinoWeight",&neutrinoWeight,"neutrinoWeight/D"); 
  RecoTree->Branch("IsItBelowStation",&IsItBelowStation,"IsItBelowStation/I");
  RecoTree->Branch("pnu",&pnu,"pnu/D");  
  RecoTree->Branch("ArrivalDirection",ArrivalDirection,"ArrivalDirection[3]"); 
  RecoTree->Branch("timeRay_AS",timeRay_AS,"timeRay_AS[2][16]/D");

  TGraph *grCor[MCH];
  TGraph *gr[MCH];
  double CorScore[MCH];
  //Start looping over the events in the file
  for(int ievt=0; ievt<nument; ievt++){
    std::cout<<"  "<<std::endl;
    std::cout<<"  "<<std::endl;
    eventTree->GetEntry(ievt);
    simTree->GetEntry(ievt);
    neutrinoWeight=weight;
    pnu= eventPtr->pnu;
    
    std::cout<<"Looking at event number "<<ievt<<" with pnu "<<pnu<<std::endl;

    //Get Event X,Y and Z from AraSim
    double event_position[3];
    event_position[0]=(eventPtr->Nu_Interaction[0].posnu.GetX());
    event_position[1]=(eventPtr->Nu_Interaction[0].posnu.GetY());
    event_position[2]=(eventPtr->Nu_Interaction[0].posnu.GetZ());
    //printf("Event Position %.2f, %.2f, %.2f \n", event_position[0], event_position[1], event_position[2]);

    //Get average coordinates for the ARA detector that was simulated with AraSim
    double *average_position = get_detector_cog(detector);
    //printf("Detector Center %.2f, %.2f, %.2f \n", average_position[0], average_position[1], average_position[2]);
  
    //initialize variables to store X,Y and Z of the vertex w.r.t the ARA detector
    TrueTxCor_XYZ[0]=event_position[0]-average_position[0];
    TrueTxCor_XYZ[1]=event_position[1]-average_position[1];
    TrueTxCor_XYZ[2]=event_position[2]-average_position[2];
    TrueTxCor_ThPhR[0]=0;
    TrueTxCor_ThPhR[1]=0;
    TrueTxCor_ThPhR[2]=0;

    //Calculate the theta, phi and R of the event w.r.t the detector center
    Interferometer::XYZtoThPhR(TrueTxCor_XYZ,TrueTxCor_ThPhR);
    TrueTxCor_ThPhR[0]=TrueTxCor_ThPhR[0]*(180.0/Interferometer::pi);
    TrueTxCor_ThPhR[1]=TrueTxCor_ThPhR[1]*(180.0/Interferometer::pi);

    //Fill in the event details
    eventNum=ievt;
    runNum=Run;

    //Define arrays which will be used later to filter out channels which were faulty or had low SNR
    double CutCh[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double FaultyCh[16]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
   
    //Variables to count the number of channels that are available after some were cut
    int NumChAvailable=0;
    int NumChAvailableV=0;
    int NumChAvailableH=0;
  
    //start looping over channels and finding peaks
    for(int ich=0; ich<MCH; ich++){
      
      //Get the Waveform from the data file for each channel
      TGraph *grdum=usefulAtriEvPtr->getGraphFromRFChan(ich);

      //Interpolate the waveforms to ensure equal spacing between the samples
      gr[ich]=FFTtools::getInterpolatedGraph(grdum,0.01);
      TGraph *grc=FFTtools::getInterpolatedGraph(grdum,0.01);
      delete grdum;

      //cout<<"total bins are "<<gr[ich]->GetN()<<endl;

      //Initiate values to set up the peak finding algorithm
      int numBinsToIntegrate=500;
      string title="ich";
      TGraph* grIntPower = makeIntegratedBinPowerGraphs(grc, numBinsToIntegrate,  title);

      vector<double> vvHitTimes; // vector of vector of hit times (first index is antenna, second index is hit time)
      vector<double> vvPeakIntPowers; // vector of vector of peak values
     
      //Find the peaks
      int numSearchPeaks=2; //only look for two peaks
      double peakSeparation=30.0;
      getAbsMaximum_N(grIntPower, numSearchPeaks, peakSeparation, vvHitTimes, vvPeakIntPowers);
      
      //Calculate the noise RMS for the waveform
      double noiserms=GetNoiseRMS(grIntPower);
      //double noiserms=1;
      
      //Fill in the peak hit times and peak SNR values
      ChHitTime[0][ich]=0;
      ChHitTime[1][ich]=0;
      ChHitTime_WF[0][ich]=vvHitTimes[0];
      ChHitTime_WF[1][ich]=vvHitTimes[1];
      ChSNR[0][ich]=vvPeakIntPowers[0]/noiserms;
      ChSNR[1][ich]=vvPeakIntPowers[1]/noiserms; 
      IgnoreCh[0][ich]=1;
      IgnoreCh[1][ich]=1;
     
      ////If the SNR of the second peak is 20 percent lower than the first one then ignore the second peak
      if(ChSNR[1][ich]/ChSNR[0][ich]<0.2){
	ChSNR[1][ich]=10;
	IgnoreCh[1][ich]=0;
      }
      
      //Ignore the second peak for now by giving it a below threshold SNR
      //ChSNR[1][ich]=10;
      
      //Filter out channels which have low SNR or are faulty
      if( (ChSNR[0][ich]<25 && ChSNR[1][ich]<25) || FaultyCh[ich]==1){
      //if(FaultyCh[ich]==1){
	CutCh[ich]=1;
      }else{
	NumChAvailable++;
	if(ich<8){
	  NumChAvailableV++;
	}else{
	  NumChAvailableH++;
	}
      }
    
      delete grc;
    }
   
    //Loop over channels and fill in the arrays for hit times and SNR that will be passed on to the interferometer
    for(int ich=0; ich<16; ich++){
    
      //Ignore channels which have been cut
      if(CutCh[ich]==-1){

    	double Dtime=0,Rtime=0;
    	double DSNR=0,RSNR=0;
	
    	//If both peaks are present
    	if(ChSNR[0][ich]>=25 && ChSNR[1][ich]>=25){
	    
    	  Dtime=ChHitTime_WF[0][ich];
    	  Rtime=ChHitTime_WF[1][ich];
    	  DSNR=ChSNR[0][ich];
    	  RSNR=ChSNR[1][ich];

    	  if(Dtime>Rtime){
    	    swap(Dtime,Rtime);
    	    swap(ChSNR[0][ich],ChSNR[1][ich]);
    	  }

    	  // if(ChSNR[0][ich]<ChSNR[1][ich]){
    	  //   swap(Dtime,Rtime);
    	  //   swap(ChSNR[0][ich],ChSNR[1][ich]);
    	  //   Rtime=0;
    	  // }
	      
    	  ChHitTime[0][ich]=Dtime;
    	  ChHitTime[1][ich]=Rtime;

    	  cout<<ich<<" two peak "<<Dtime<<" "<<Rtime<<" "<<ChSNR[0][ich]<<" "<<ChSNR[1][ich]<<endl;
    	}else{
    	  //If only one peak is present

    	  //Dtime=ChHitTime[0][ich];
    	  Dtime=ChHitTime_WF[0][ich];

    	  ChHitTime[0][ich]=Dtime;
    	  ChHitTime[1][ich]=0;
    	  ChSNR[1][ich]=0;
    	  cout<<ich<<" one peak "<<Dtime<<" "<<ChSNR[0][ich]<<endl;
    	}

      }//cut channel if condition

    }////channel loop
  
    Double_t NumSinglePeak=0;
    Double_t NumTotalChannels=0;
    vector <double> hittimes;
    //Loop over channels and fill in the ignore channel array which will be used by the interferometer to see which channels need to be ignored
    for(int ich=0;ich<TotalAntennasRx;ich++){
      for(int iray=0;iray<2;iray++){
	
    	if(fabs(ChHitTime[iray][ich])<1e-5 || CutCh[ich]==1){
    	  IgnoreCh[iray][ich]=0;
    	}
      
      }
      if(IgnoreCh[0][ich]==1){
	hittimes.push_back(ChHitTime[0][ich]);
    	NumTotalChannels++;
      }
      if((IgnoreCh[0][ich]==1 && IgnoreCh[1][ich]==0)){
    	NumSinglePeak++;
      }
    } 

    //If in all the selected channels 85% have single peak and 15% have double peaks then ignore the second peak in the channels having double peaks. This is probably because those second peaks are noise peaks.
    Double_t FractionSinglePeak=NumSinglePeak/NumTotalChannels;
    if(FractionSinglePeak>=0.85){
      for(int ich=0;ich<TotalAntennasRx;ich++){
    	if(IgnoreCh[1][ich]==1){
    	  IgnoreCh[1][ich]=0;
    	  if(ChSNR[1][ich]>ChSNR[0][ich]){
    	    swap(ChHitTime[0][ich], ChHitTime[1][ich]);
    	    swap(ChSNR[0][ich], ChSNR[1][ich]);
    	  }
    	}
      
      }
    }
    
    //Print out the peak hit times for all channels and see what peaks have been passed through and what peaks have been ignored. If IgnoreCh is 1 then peak has passed and if IgnoreCh  0 then the peak has been ignored.
    // for(int ich=0;ich<16;ich++){
    //   cout<<ich<<" "<<IgnoreCh[0][ich]<<" "<<ChHitTime[0][ich]<<" "<<IgnoreCh[1][ich]<<" "<<ChHitTime[1][ich]<<endl;
    // }

    IsItBelowStation=-1;// this variable is used to check if station was hit from above or below. Just dummy fill it for now.

    int IgnoreChB[2][16];
  
    for(int ich=0;ich<16;ich++){
      IgnoreChB[0][ich]=1;
      IgnoreChB[1][ich]=1;
      for(int iray=0;iray<2;iray++){
    	timeRay_AS[iray][ich]=-1000;
      }
    }

    //start looping over channels and finding peaks
    for(int ist=0; ist<4; ist++){
      for(int iant=0; iant<4; iant++){
	int rfchan=getRFChannel(ist, iant);
	if(reportPtr->stations[0].strings[ist].antennas[iant].arrival_time.size()==2){
	  timeRay_AS[0][rfchan]=reportPtr->stations[0].strings[ist].antennas[iant].arrival_time[0]*pow(10,9);
	  timeRay_AS[1][rfchan]=reportPtr->stations[0].strings[ist].antennas[iant].arrival_time[1]*pow(10,9);
	  //cout<<ist<<" "<<iant<<" "<<rfchan<<endl;
	  if(timeRay_AS[0][rfchan]>timeRay_AS[1][rfchan]){
	    swap(timeRay_AS[0][rfchan],timeRay_AS[1][rfchan]);
	  }
	}

	if(reportPtr->stations[0].strings[ist].antennas[iant].arrival_time.size()==1){
	  timeRay_AS[0][rfchan]=reportPtr->stations[0].strings[ist].antennas[iant].arrival_time[0]*pow(10,9);
	  timeRay_AS[1][rfchan]=-1000;
	}

	if(reportPtr->stations[0].strings[ist].antennas[iant].arrival_time.size()==0){
	  timeRay_AS[0][rfchan]=-1000;
	  timeRay_AS[1][rfchan]=-1000;
	}

      }
    }    

    vector <double> timeRayv[2];
    vector <int> IgnoreChBv[2];
    vector <int> IgnoreChv[2];
    vector <double> ChSNRv[2];
    timeRayv[0].resize(TotalAntennasRx);
    timeRayv[1].resize(TotalAntennasRx);
   
    IgnoreChBv[0].resize(TotalAntennasRx);
    IgnoreChBv[1].resize(TotalAntennasRx);
    IgnoreChv[0].resize(TotalAntennasRx);
    IgnoreChv[1].resize(TotalAntennasRx);
    ChSNRv[0].resize(TotalAntennasRx);
    ChSNRv[1].resize(TotalAntennasRx);

    for(int ich=0;ich<MCH;ich++){
      IgnoreChBv[0][ich]=1;
      IgnoreChBv[1][ich]=1;
      IgnoreChv[0][ich]=1;
      IgnoreChv[1][ich]=1;
      ChSNRv[0][ich]=1;
      ChSNRv[1][ich]=1;
      timeRayv[0][ich]=0;
      timeRayv[1][ich]=0;
    }

    vector <double> ChHitTimev[2];
    ChHitTimev[0].resize(TotalAntennasRx);
    ChHitTimev[1].resize(TotalAntennasRx);

    for(int ich=0;ich<TotalAntennasRx;ich++){
      for(int iray=0;iray<2;iray++){
	ChHitTimev[iray][ich]=ChHitTime[iray][ich];
	IgnoreChv[iray][ich]=IgnoreCh[iray][ich];
	ChSNRv[iray][ich]=ChSNR[iray][ich];
      }
      // if(fabs(ChHitTime[0][ich] - ChHitTime[1][ich])>1000){
      // 	IgnoreCh[1][ich]=0;
      // 	IgnoreChv[1][ich]=0;
      // }
    }

    //Only do reconstruction if more than 4 channels are available for reconstruction
    bool CheckRRays=false;
    if(NumChAvailable>=4){
      
      //Check if the station was hit from above or below
      IsItBelowStation=Interferometer::IsItAboveOrBelow(ChHitTimev,IgnoreChv) ;
    
      for(int ich=0;ich<16;ich++){
      	cout<<"check again "<<ich<<" "<<IgnoreChv[0][ich]<<" "<<ChHitTimev[0][ich]<<" "<<ChSNRv[0][ich]<<" "<<IgnoreChv[1][ich]<<" "<<ChHitTimev[1][ich]<<" "<<ChSNRv[1][ich]<<endl;
	if(IgnoreCh[1][ich]!=0){
	  CheckRRays=true;
	  IsItBelowStation=1;
	}
      }

      // double ThPhRa[3]={107.039 ,5.65878,1742.06};
      // double X2a=Interferometer::GetChiSquaredThPhR(ThPhRa, ChHitTime, IgnoreCh, ChSNR, IsItBelowStation,0,50);
      // cout<<"chi sqaured values are "<<X2a<<" "<<0<<endl;

      Interferometer::GetRecieveAngle(ChHitTimev, IgnoreChv, ChSNRv, 100, ArrivalDirection);

      //print out the true vertex values of theta, phi and R 
      cout<<"True values are: |  Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
      //cout<<"True values are: |  Z_true="<<TrueTxCor_XYZ[0]<<" ,Y_true="<<TrueTxCor_XYZ[1]<<" ,Z_true="<<TrueTxCor_XYZ[2]-180<<endl;
      
      //Set up variables for interferometry
      double GuessResultCor[3][4]; 

      //Start counting time it takes to guess the initial condition
      auto t1 = std::chrono::high_resolution_clock::now();

      ////Get a guess direction
      Interferometer::GetApproximateMinThPhR(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,10);      

      //Get a guess distance
      Interferometer::GetApproximateDistance(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,50);
      
      //If events have alot of channels with noise in them it can often causes the guess value of the distance to be around 50 m which is really close to the station and highly unlikely. So in this case I go through the channels again and remove any channels which have an SNR lower than 36 in power. Then I count the number of remaining channels left to see if we have channels left to do a reconstruction
      if(GuessResultCor[0][2]<60 && IsItBelowStation==1){
      	NumChAvailable=0;
      	for(int ich=0;ich<TotalAntennasRx;ich++){
	     
      	  if(ChSNRv[0][ich]<36 && IgnoreChv[0][ich]==1){
      	    IgnoreChv[0][ich]=0;
      	    ChHitTimev[0][ich]=0;	
      	  }
      	  if(IgnoreChv[0][ich]==1){
      	    NumChAvailable++;
      	  }	   
      
      	}

      	if(NumChAvailable>=4){
      	  Interferometer::GetApproximateDistance(GuessResultCor,ChHitTimev,IgnoreChv,ChSNRv,IsItBelowStation,50);
      	}
      }

      //Calculate the time it took to guess the initial condition
      auto t2 = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
      DurationInitialCondition=duration/1000;      

      if(NumChAvailable>=4){
       
      	//Convert initial guess to radians
      	InitialTxCor_ThPhR[0]=GuessResultCor[0][0]*(Interferometer::pi/180);
      	InitialTxCor_ThPhR[1]=GuessResultCor[0][1]*(Interferometer::pi/180);
      	InitialTxCor_ThPhR[2]=GuessResultCor[0][2];  

      	//Set the radial width around guess distance value in meters    
      	double MinimizerRadialWidth=500;

      	//do the interferometry
      	Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTimev, IgnoreChv, ChSNRv, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   
      	cout<<"1st try Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]*(180.0/Interferometer::pi)<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]*(180.0/Interferometer::pi)<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;

      	InitialTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
      	InitialTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180);
      	InitialTxCor_ThPhR[2]=FinalTxCor_ThPhR[2];  

      	Interferometer::DoInterferometery(InitialTxCor_ThPhR, FinalTxCor_ThPhR, ChHitTimev, IgnoreChv, ChSNRv, FinalMinValue, DurationReconstruction, Iterations,MinimizerRadialWidth, IsItBelowStation,50);   
    
      	//Calcualte the total duration
      	DurationTotal=DurationInitialCondition+DurationReconstruction;
    
      	//Convert Initial guess direction from theta,phi,r to x,y,z
      	InitialTxCor_XYZ[0]=0;
      	InitialTxCor_XYZ[1]=0;
      	InitialTxCor_XYZ[2]=0;
      	Interferometer::ThPhRtoXYZ(InitialTxCor_ThPhR,InitialTxCor_XYZ);
      	InitialTxCor_ThPhR[0]=InitialTxCor_ThPhR[0]*(180./Interferometer::pi);
      	InitialTxCor_ThPhR[1]=InitialTxCor_ThPhR[1]*(180./Interferometer::pi); 

      	//Convert final reconstructed direction from theta,phi,r to x,y,z
      	FinalTxCor_XYZ[0]=0;
      	FinalTxCor_XYZ[1]=0;
      	FinalTxCor_XYZ[2]=0;
      	FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(Interferometer::pi/180);
      	FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(Interferometer::pi/180); 
      	Interferometer::ThPhRtoXYZ(FinalTxCor_ThPhR, FinalTxCor_XYZ);
      	FinalTxCor_ThPhR[0]=FinalTxCor_ThPhR[0]*(180./Interferometer::pi);
      	FinalTxCor_ThPhR[1]=FinalTxCor_ThPhR[1]*(180./Interferometer::pi);
  
      	//cout<<"Final Reco Results are: |  X_initial="<<InitialTxCor_XYZ[0]<<" ,Y_initial="<<InitialTxCor_XYZ[1]<<" ,Z_initial="<<InitialTxCor_XYZ[2]<<" | X_reco="<<FinalTxCor_XYZ[0]<<" ,Y_reco="<<FinalTxCor_XYZ[1]<<" ,Z_reco="<<FinalTxCor_XYZ[2]<<endl;
      	cout<<"Final Reco Results are: |  Th_initial="<<InitialTxCor_ThPhR[0]<<" ,Ph_initial="<<InitialTxCor_ThPhR[1]<<" ,R_initial="<<InitialTxCor_ThPhR[2]<<" | Th_reco="<<FinalTxCor_ThPhR[0]<<" ,Ph_reco="<<FinalTxCor_ThPhR[1]<<" ,R_reco="<<FinalTxCor_ThPhR[2]<<endl;
      	cout<<"True values are: |  Th_true="<<TrueTxCor_ThPhR[0]<<" ,Ph_true="<<TrueTxCor_ThPhR[1]<<" ,R_true="<<TrueTxCor_ThPhR[2]<<endl;
      	cout<<"Fn Min Value:"<<FinalMinValue<<", Total Minimizer Iterations: "<<Iterations<<", Total Reco Duration (ms): "<<DurationTotal<<endl;
    }else{

    	//If no reconstruction was performed then fill the final arrays with zeroes
    	ArrivalDirection[0]=0;
    	ArrivalDirection[1]=0;
    	ArrivalDirection[2]=0;
	  
    	FinalTxCor_XYZ[0]=0;
    	FinalTxCor_XYZ[1]=0;
    	FinalTxCor_XYZ[2]=0;
	  
    	FinalTxCor_ThPhR[0]=0;
    	FinalTxCor_ThPhR[1]=0;
    	FinalTxCor_ThPhR[2]=0;
	  
    	cout<<" Number of hit channels is less than 3 so no reconstruction was performed!!!"<<endl;
      }   
	

    }else{

      //If no reconstruction was performed then fill the final arrays with zeroes
      ArrivalDirection[0]=0;
      ArrivalDirection[1]=0;
      ArrivalDirection[2]=0;

      FinalTxCor_XYZ[0]=0;
      FinalTxCor_XYZ[1]=0;
      FinalTxCor_XYZ[2]=0;

      FinalTxCor_ThPhR[0]=0;
      FinalTxCor_ThPhR[1]=0;
      FinalTxCor_ThPhR[2]=0;

      cout<<" Number of hit channels is less than 3 so no reconstruction was performed!!!"<<endl;
    }

    //Fill the tree and write all variables to it
    OutputFile->cd();
    RecoTree->Fill();
    RecoTree->Write();

  }  

  delete InputFile;

  OutputFile->Write();
  OutputFile->Close(); 
  
}
