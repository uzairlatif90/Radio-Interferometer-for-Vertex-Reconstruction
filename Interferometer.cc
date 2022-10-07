#include "Interferometer.hh"

void Interferometer::XYZtoThPhR(Double_t XYZ[3],Double_t ThPhR[3]){
  ThPhR[0]=atan2(sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]),XYZ[2]);
  ThPhR[1]=atan2(XYZ[1],XYZ[0]);
  ThPhR[2]=sqrt(XYZ[0]*XYZ[0] + XYZ[1]*XYZ[1] + XYZ[2]*XYZ[2]);
}

void Interferometer::ThPhRtoXYZ(Double_t ThPhR[3],Double_t XYZ[3]){
  XYZ[0]=ThPhR[2]*sin(ThPhR[0])*cos(ThPhR[1]);
  XYZ[1]=ThPhR[2]*sin(ThPhR[0])*sin(ThPhR[1]);
  XYZ[2]=ThPhR[2]*cos(ThPhR[0]);
}

void Interferometer::GenerateChHitTimeAndCheckHits(double TxCor[3],vector<double> (&timeRay)[2],vector<int> (&IgnoreCh)[2]){
  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);

  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]+AvgAntennaCoordRx[2]};
  
  int DHits=0,RHits=0;
  double timeD, timeR, timeRa[2];
  double RangD, RangR, RangRa[2];
  double LangD, LangR, LangRa[2];

  vector <double> RangRay[2];  
  vector <double> LangRay[2];  
  RangRay[0].resize(TotalAntennasRx);
  LangRay[0].resize(TotalAntennasRx);
  RangRay[1].resize(TotalAntennasRx);
  LangRay[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[0][iRx]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[1][iRx]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[2][iRx]+AvgAntennaCoordRx[2];
    double TxDepth=AntennaCoordTx[2];
    
    //cout<<"rt parameters are "<<0<<" "<<TxDepth<<" "<<Distance<<" "<<RxDepth<<endl;
    //cout<<"A B C "<<IceRayTracing::A_ice<<" "<<IceRayTracing::B_ice<<" "<<IceRayTracing::C_ice<<endl;

    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);
   
    // cout<<"*******For the Direct Ray********"<<endl;
    // cout<<"Launch Angle: "<<RTresults[0]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<RTresults[8]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<RTresults[4]*pow(10,9)<<" ns"<<endl;
    // cout<<"*******For the Reflected Ray********"<<endl;
    // cout<<"Launch Angle: "<<RTresults[1]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<RTresults[9]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<RTresults[5]*pow(10,9)<<" ns"<<endl;   
    // cout<<"Incident Angle in Ice on the Surface: "<<RTresults[18]<<" deg"<<endl;
    // cout<<"*******For the Refracted[1] Ray********"<<endl;
    // cout<<"Launch Angle: "<<RTresults[2]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<RTresults[10]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<RTresults[6]*pow(10,9)<<" ns"<<endl;
    // cout<<"*******For the Refracted[2] Ray********"<<endl;
    // cout<<"Launch Angle: "<<RTresults[3]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<RTresults[11]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<RTresults[7]*pow(10,9)<<" ns"<<endl;
    
    timeD=RTresults[4]*pow(10,9);
    timeR=RTresults[5]*pow(10,9);
    timeRa[0]=RTresults[6]*pow(10,9);
    timeRa[1]=RTresults[7]*pow(10,9);

    RangD=RTresults[8];
    RangR=RTresults[9];
    RangRa[0]=RTresults[10];
    RangRa[1]=RTresults[11];
    
    timeRay[0][iRx]=timeD;
    timeRay[1][iRx]=timeR;

    RangRay[0][iRx]=RangD;
    RangRay[1][iRx]=RangR;

    if(RangD!=-1000){
      timeRay[0][iRx]=timeD;
      RangRay[0][iRx]=RangD;
      LangRay[0][iRx]=LangD;
    }

    if(RangR!=-1000){
      timeRay[1][iRx]=timeR;
      RangRay[1][iRx]=RangR;
      LangRay[1][iRx]=LangR;
    }
    
    if(RangRa[0]!=-1000 && RangD!=-1000){
      timeRay[0][iRx]=timeD;
      RangRay[0][iRx]=RangD;
      LangRay[0][iRx]=LangD;
      timeRay[1][iRx]=timeRa[0];
      RangRay[1][iRx]=RangRa[0];
      LangRay[1][iRx]=LangRa[0];
    }
      
    if(RangRa[0]!=-1000 && RangR!=-1000){
      timeRay[1][iRx]=timeR;
      RangRay[1][iRx]=RangR;
      LangRay[1][iRx]=LangR;
      timeRay[0][iRx]=timeRa[0];
      RangRay[0][iRx]=RangRa[0];
      LangRay[0][iRx]=LangRa[0];
    }

    if(RangRa[1]!=-1000 && RangD!=-1000){
      timeRay[0][iRx]=timeD;
      RangRay[0][iRx]=RangD;
      LangRay[0][iRx]=LangD;
      timeRay[1][iRx]=timeRa[1];
      RangRay[1][iRx]=RangRa[1];
      LangRay[1][iRx]=LangRa[1];
    }
      
    if(RangRa[1]!=-1000 && RangR!=-1000){
      timeRay[1][iRx]=timeR;
      RangRay[1][iRx]=RangR;
      LangRay[1][iRx]=LangR;
      timeRay[0][iRx]=timeRa[1];
      RangRay[0][iRx]=RangRa[1];
      LangRay[0][iRx]=LangRa[1];
    }

    if(RangRa[1]!=-1000 && RangRa[0]!=-1000){
      timeRay[1][iRx]=timeRa[1];
      RangRay[1][iRx]=RangRa[1];
      LangRay[1][iRx]=LangRa[1];
      timeRay[0][iRx]=timeRa[0];
      RangRay[0][iRx]=RangRa[0];
      LangRay[0][iRx]=LangRa[0];
    }

    if(RangRay[1][iRx]==-1000 && RangRay[0][iRx]==-1000 && RangRa[0]!=-1000){
      timeRay[0][iRx]=timeRa[0];
      RangRay[0][iRx]=RangRa[0];
      LangRay[0][iRx]=LangRa[0];
    }

    if(RangRay[1][iRx]==-1000 && RangRay[0][iRx]==-1000 && RangRa[1]!=-1000){
      timeRay[1][iRx]=timeRa[1];
      RangRay[1][iRx]=RangRa[1];
      LangRay[1][iRx]=LangRa[1];
    }
    
    if(RangRay[0][iRx]==-1000){
      IgnoreCh[0][iRx]=0;
    }

    if(RangRay[1][iRx]==-1000){
      IgnoreCh[1][iRx]=0;
    }

    if(timeRay[0][iRx]>timeRay[1][iRx] && RangRay[0][iRx]!=-1000 && RangRay[1][iRx]!=-1000){
      swap(timeRay[0][iRx],timeRay[1][iRx]);
    }
    
    if(RangRay[0][iRx]!=-1000){
      DHits++;
    }

    if(RangRay[1][iRx]!=-1000){
      RHits++;
    }
    delete [] RTresults;  
  }  
  //cout<<"hits are "<<DHits<<" "<<RHits<<endl; 
}

void Interferometer::GenerateChHitTimeAndCheckHits_Cnz(double TxCor[3],vector<double> (&timeRay)[2],vector<int> (&IgnoreCh)[2]){

  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  
  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]+AvgAntennaCoordRx[2]};
  
  int DHits=0,RHits=0;

  double timeD, timeR, timeRa;
  double RangD, RangR, RangRa;
  vector <double> RangRay[2];  
  RangRay[0].resize(TotalAntennasRx);
  RangRay[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[0][iRx]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[1][iRx]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[2][iRx]+AvgAntennaCoordRx[2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing_Cnz(0, TxDepth, Distance,RxDepth,(IceRayTracing::Getnz(TxDepth)+IceRayTracing::Getnz(RxDepth))/2);

    timeD=RTresults[2]*pow(10,9);
    timeR=RTresults[3]*pow(10,9);
   
    RangD=RTresults[4];
    RangR=RTresults[5];
    
    timeRay[0][iRx]=timeD;
    timeRay[1][iRx]=timeR;

    RangRay[0][iRx]=RangD;
    RangRay[1][iRx]=RangR;
    
    DHits++;
    RHits++;
    
    delete [] RTresults;  
  }  
  //cout<<"hits are "<<DHits<<" "<<RHits<<endl; 
}


void Interferometer::GenerateChHitTimeAndCheckHits_Air(double TxCor[3],vector<double> (&timeRay)[2],vector<int> (&IgnoreCh)[2]){

  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  
  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]+AvgAntennaCoordRx[2]};
  
  int DHits=0;
  double timeD;
  double RangD;
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[0][iRx]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[1][iRx]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[2][iRx]+AvgAntennaCoordRx[2];
    double TxHeight=AntennaCoordTx[2];   
    //cout<<"RT par are "<<RxDepth<<" "<<Distance<<" "<<TxHeight<<endl;
    double *RTresults=IceRayTracing::GetDirectRayPar_Air(RxDepth, Distance, TxHeight);
    timeD=RTresults[2]*pow(10,9);
    RangD=RTresults[0];
    timeRay[0][iRx]=timeD;

    if(RangD==-1000){
      IgnoreCh[0][iRx]=0;
    }

    IgnoreCh[1][iRx]=0;
    
    DHits++;   
    
    delete [] RTresults;  
  }
}

bool Interferometer::CheckTrigger(vector<int> (&IgnoreCh)[2]){

  int count[2][2]={{0,0},{0,0}};
  
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(IgnoreCh[iray][iRx]!=0 && iRx<8){
	count[iray][0]++;
      }
      if(IgnoreCh[iray][iRx]!=0 && iRx>7){
	count[iray][1]++;
      }
    }
  }
  
  bool RayTrigger[2]={false,false};

  for(int iray=0;iray<2;iray++){
    if(count[iray][0]+count[iray][1]>=4){
      if(count[iray][0]>=4){
	RayTrigger[iray]=true;
      } 
      if(count[iray][1]>=4){
	RayTrigger[iray]=true;
      }
    }
  }

  //cout<<"counts are "<<count[0][0]<<" "<<count[0][1]<<" "<<count[1][0]<<" "<<count[1][1]<<endl;

  bool DidStationTrigger=false;

  if(RayTrigger[0]==true || RayTrigger[1]==true){
    DidStationTrigger=true;
  }

  if(RayTrigger[0]==false && RayTrigger[1]==false){
    DidStationTrigger=false;
  }

  return DidStationTrigger;
}

void Interferometer::ReadChHitTimeFromData(const char * filename,vector<double> (&ChHitTime)[2]){
  
  ChHitTime[0].resize(TotalAntennasRx);
  ChHitTime[1].resize(TotalAntennasRx);
  
  ////Open the file
  std::ifstream ain(filename);
  int n1=0;////variable for counting total number of data points
  std::string line;
  double dummy[3]={0,0,0};////temporary variable for storing data values from the file  

  //Check if file is open and store data
  if(ain.is_open()){
    while (true){
      ain>>dummy[0]>>dummy[1]>>dummy[2];
      if( ain.eof() ) break;
      ChHitTime[0][n1]= dummy[1];
      ChHitTime[1][n1]= dummy[2];
      n1++;
    }
  }
  ain.close();
}

void Interferometer::FindFirstHitAndNormalizeHitTime(vector<double> (&ChHitTime)[2],vector<int> (&IgnoreCh)[2], int NormalizeRay){

  double FirstHitTime[2]={0,0};
  int FirstHitCh[2]={0,0};

  double LastHitTime[2]={0,0};
  int LastHitCh[2]={0,0};
  vector <double> ChHitTimeB[2];
  vector <int> ChHitNum[2];
  
  for(int iray=0;iray<2;iray++){ 
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(IgnoreCh[iray][iRx]==1){
	ChHitTimeB[iray].push_back(ChHitTime[iray][iRx]);
	ChHitNum[iray].push_back(iRx);
      }
    }
  }
  
  for(int iray=0;iray<2;iray++){
    if(ChHitTimeB[iray].size()!=0){
      FirstHitTime[iray]=TMath::MinElement(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      FirstHitCh[iray]=TMath::LocMin(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      
      LastHitTime[iray]=TMath::MaxElement(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
      LastHitCh[iray]=TMath::LocMax(ChHitTimeB[iray].size(),ChHitTimeB[iray].data());
    }
  }
  
  for(int iray=0;iray<2;iray++){
    if(ChHitTimeB[iray].size()!=0){
      for(int iRx=0;iRx<TotalAntennasRx;iRx++){
	if(IgnoreCh[iray][iRx]==1){
	  ChHitTime[iray][iRx]=ChHitTime[iray][iRx]-FirstHitTime[NormalizeRay];
	}else{
	  ChHitTime[iray][iRx]=0;
	}
      }
    }else{
      for(int iRx=0;iRx<TotalAntennasRx;iRx++){
	ChHitTime[iray][iRx]=0;
      }
    }
  }

  // for(int iRx=0;iRx<TotalAntennasRx;iRx++){
  //   //ChHitTime[iray][iRx]=0;
  //   cout<<iRx<<" "<<ChHitTime[0][iRx]<<" "<<ChHitTime[0][iRx]<<endl;
  // }
}

void Interferometer::AddGaussianJitterToHitTimes(double JitterNumber,vector<double> (&ChHitTime)[2]){

  ////Introduce jitter +-JitterNumber
  TRandom3 *RandNumIni = new TRandom3(0);
  for (int iRx=0; iRx<TotalAntennasRx; iRx++){
    for(int iray=0;iray<2;iray++){ 
      double RandNum = (RandNumIni->Rndm(iRx)*2*JitterNumber)-JitterNumber;
      //double RandNum = (RandNumIni->Gaus(0,0.4));
      if(RandNum>=JitterNumber){
      	RandNum=JitterNumber;
      }
      if(RandNum<=-JitterNumber){
      	RandNum=-JitterNumber;
      }
      //cout<<"rand num is "<<RandNum<<endl;
      ChHitTime[iray][iRx]=ChHitTime[iray][iRx]+RandNum;
    }
  }
  delete RandNumIni;
}

int Interferometer::IsItAboveOrBelow(vector<double> (&ChHitTime)[2],vector<int> (&IgnoreCh)[2]){

  int IsItBelowStation=1;
  double XYZ[3]={0,0,0};
  double min=1e9;
  int min_Ch=0;
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(IgnoreCh[0][iRx]==1){
      if(ChHitTime[0][iRx]<min){
	min=ChHitTime[0][iRx];
	min_Ch=iRx;
      }
    }
  }
  
  // XYZ[0]=AntennaCoordRx[min_Ch][0]-AvgAntennaCoordRx[0];
  // XYZ[1]=AntennaCoordRx[min_Ch][1]-AvgAntennaCoordRx[1];
  // XYZ[2]=AntennaCoordRx[min_Ch][2]-AvgAntennaCoordRx[2];

  XYZ[0]=AntennaCoordRx[0][min_Ch];
  XYZ[1]=AntennaCoordRx[1][min_Ch];
  XYZ[2]=AntennaCoordRx[2][min_Ch];

  if(XYZ[2]>0){
    cout<<"Ch."<<min_Ch<<" was hit first and the arrival direction is from above the station"<<endl;
    IsItBelowStation=0;
  }
  if(XYZ[2]<0){
    cout<<"Ch."<<min_Ch<<" was hit first and the arrival direction is from below the station in ice"<<endl;
    IsItBelowStation=1;
  }
   
  return IsItBelowStation;
}

double Interferometer::Minimizer_f1D(double m, void *params){

  double *p = (double *)params;
  
  double ExpectedX=p[6*TotalAntennasRx];
  double ExpectedY=p[6*TotalAntennasRx+1];
  double ExpectedZ=p[6*TotalAntennasRx+2];

  double ExpectedTh=p[6*TotalAntennasRx+3];
  double ExpectedPh=p[6*TotalAntennasRx+4];
  double ExpectedR=p[6*TotalAntennasRx+5];
  
  double theta, phi, r;
  theta = ExpectedTh*(Interferometer::pi/180.0);
  phi = ExpectedPh*(Interferometer::pi/180.0);
  r = m;
  
  Double_t ThPhR[3]={theta,phi,r};
  Double_t XYZ[3]={0,0,0};
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
  
  double output=0;

  double AntennaCoordTx[3]={0,0,0};
  for(int ixyz=0;ixyz<3;ixyz++){
    AntennaCoordTx[ixyz]=XYZ[ixyz];
  }
    
  vector<double> timeRay[2];
  vector<int> IgnoreCh[2];

  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }    

  if(XYZ[2]+AvgAntennaCoordRx[2]<0 || p[6*TotalAntennasRx+12]==0){
    //cout<<"working for ice "<<ExpectedTh<<" "<<ExpectedPh<<" "<<p[6*TotalAntennasRx+10]<<" "<<XYZ[2]+AvgAntennaCoordRx[2]<<" "<<p[6*TotalAntennasRx+12]<<endl;;
    Interferometer::GenerateChHitTimeAndCheckHits(AntennaCoordTx,timeRay,IgnoreCh);  
  }else{
    //cout<<"working for air "<<ExpectedTh<<" "<<ExpectedPh<<" "<<p[6*TotalAntennasRx+10]<<" "<<XYZ[2]+AvgAntennaCoordRx[2]<<" "<<p[6*TotalAntennasRx+12]<<endl;;
    Interferometer::GenerateChHitTimeAndCheckHits_Air(AntennaCoordTx,timeRay,IgnoreCh);  
  }
  
  Double_t SumSNRD=0, SumSNRR=0;
  int chanD[2]={0,0};
  int chanR[2]={0,0};
  int chanDsame=0,chanRsame=0,chanRsameB=0;
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(p[2*TotalAntennasRx +iRx]!=0){
      chanD[0]++;
    }
    if(p[3*TotalAntennasRx +iRx]!=0){
      chanR[0]++;
    }
    if(IgnoreCh[0][iRx]!=0){
      chanD[1]++;
    }
    if(IgnoreCh[1][iRx]!=0){
      chanR[1]++;
    }
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
      SumSNRD+=p[4*TotalAntennasRx +iRx];
      chanDsame++;
    }
    if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      SumSNRR+=p[5*TotalAntennasRx +iRx];
      chanRsame++;
    }

    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      chanRsameB++;
    }

  }
 
  double deltaT=0;
  double output1=0;  
  double output2=0;
  
  //if((p[6*TotalAntennasRx+13]==1 && chanDsame>=chanD[0]) || ( p[6*TotalAntennasRx+13]==0 && chanRsameB>=chanD[0])){
    if(chanDsame>=chanD[0]){
      //if(p[6*TotalAntennasRx+13]==1){

      double loop1_D=0,loop1_R=0;   
      double loop2_D=0,loop2_R=0;
      double sum1_D=0,sum1_R=0;  

      for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
     
	loop2_D=0,loop2_R=0;
	double sum_D=0,sum_R=0;
	double sum_W_D=0,sum_W_R=0;

	for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 
	  if(p[2*TotalAntennasRx +iRx]!=0 && p[2*TotalAntennasRx +iRxB]!=0 && IgnoreCh[0][iRx]!=0 && IgnoreCh[0][iRxB]!=0){     

	    if(sum_D==0){
	      loop1_D+=p[4*TotalAntennasRx +iRx];
	    }
	    
	    deltaT=(timeRay[0][iRx] - timeRay[0][iRxB] - (p[0+iRx] - p[0+iRxB]));
	    sum_D+=pow(deltaT,2)*p[4*TotalAntennasRx +iRxB];
	    sum_W_D+=p[4*TotalAntennasRx +iRxB];
	    loop2_D++;

	  }
	}

	if(loop2_D>0){
	  sum_D=sum_D/sum_W_D;
	  if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
	    sum1_D+=sum_D*p[4*TotalAntennasRx +iRx];
	  }
	}

	loop2_R=0;
	sum_R=0;
	sum_W_R=0;

	for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 
	  if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0  && p[3*TotalAntennasRx +iRxB]!=0 && IgnoreCh[1][iRxB]!=0){
	    if(sum_R==0){
	      loop1_D+=p[5*TotalAntennasRx +iRx];
	      sum_R=1;
	    }

	    deltaT=(timeRay[1][iRx] - timeRay[1][iRxB] - (p[1*TotalAntennasRx+iRx] - p[1*TotalAntennasRx+iRxB]));
	    sum_D+=pow(deltaT,2)*p[5*TotalAntennasRx +iRxB];
	    sum_W_D+=p[5*TotalAntennasRx +iRxB];
	    loop2_R++;
	  }
	}

	if(loop2_R>0){
	  sum_R=sum_D/sum_W_D;
	  if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
	    sum1_D+=sum_R*p[5*TotalAntennasRx +iRx];
	  }
	}

      }
      
      output1+=sum1_D/loop1_D;
      //}  

    //if(p[6*TotalAntennasRx+13]==0){

    //   double loop1=0;
    //   double loop2=0;
    
    //   double sum1=0;  
    //   for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
    // 	loop2=0;
    // 	double sum=0;
    // 	double sum_W=0;
    // 	for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 

    // 	  if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0 && p[2*TotalAntennasRx +iRxB]!=0 && IgnoreCh[1][iRxB]!=0){

    // 	    if(sum==0){
    // 	      loop1+=p[4*TotalAntennasRx +iRx];
    // 	    }

    // 	    deltaT=(timeRay[1][iRx] - timeRay[1][iRxB]  - (p[0+iRx] - p[0+iRxB]));
    // 	    sum+=pow(deltaT,2)*p[4*TotalAntennasRx +iRxB];
    // 	    sum_W+=p[4*TotalAntennasRx +iRxB];	  
    // 	    loop2++;

    // 	  }
    // 	}
      
    // 	if(loop2>0){
    // 	  sum=sum/sum_W;
    // 	  if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
    // 	    sum1+=sum*p[4*TotalAntennasRx +iRx];
    // 	  }
    // 	}
    //   }
      
    //   output2+=sum1/loop1;

    // }
    
    //if(p[6*TotalAntennasRx+13]==1){
      output=output1;
      //}
 
    // if(p[6*TotalAntennasRx+13]==0){
    //   output=output2;
    // }
    

    }else{
    output=1e9+r;
  }

  if(std::isnan(output)){
    output=1e9+r;
  }
  
  if(ExpectedTh>179.9 || ExpectedTh<0.1){
    output=1e9+output;
  }

  if(XYZ[2]+AvgAntennaCoordRx[2]>0 && output<1e9){
    output=1e9+output;
  }

  return output;
}


double Interferometer::Minimizer_f(const gsl_vector *v, void *params){

  double *p = (double *)params;
  
  double ExpectedX=p[6*TotalAntennasRx];
  double ExpectedY=p[6*TotalAntennasRx+1];
  double ExpectedZ=p[6*TotalAntennasRx+2];

  double ExpectedTh=p[6*TotalAntennasRx+3];
  double ExpectedPh=p[6*TotalAntennasRx+4];
  double ExpectedR=p[6*TotalAntennasRx+5];
  
  double theta, phi, r;
  theta = gsl_vector_get(v, 0)*(Interferometer::pi/180.0);
  phi = gsl_vector_get(v, 1)*(Interferometer::pi/180.0);
  r = p[6*TotalAntennasRx+10];
  
  Double_t ThPhR[3]={theta,phi,r};
  Double_t XYZ[3]={0,0,0};
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
  
  double output=0;

  double AntennaCoordTx[3]={0,0,0};
  for(int ixyz=0;ixyz<3;ixyz++){
    AntennaCoordTx[ixyz]=XYZ[ixyz];
  }
    
  vector<double> timeRay[2];
  vector<int> IgnoreCh[2];

  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }    

  if(XYZ[2]+AvgAntennaCoordRx[2]<0 || p[6*TotalAntennasRx+12]==0){
    //cout<<"working for ice "<<gsl_vector_get(v, 0)<<" "<<gsl_vector_get(v, 1)<<" "<<p[6*TotalAntennasRx+10]<<" "<<XYZ[2]+AvgAntennaCoordRx[2]<<" "<<p[6*TotalAntennasRx+12]<<endl;;
    Interferometer::GenerateChHitTimeAndCheckHits(AntennaCoordTx,timeRay,IgnoreCh);  
  }else{
    //cout<<"working for air "<<gsl_vector_get(v, 0)<<" "<<gsl_vector_get(v, 1)<<" "<<p[6*TotalAntennasRx+10]<<" "<<XYZ[2]+AvgAntennaCoordRx[2]<<" "<<p[6*TotalAntennasRx+12]<<endl;;
    Interferometer::GenerateChHitTimeAndCheckHits_Air(AntennaCoordTx,timeRay,IgnoreCh);  
  }
  
  Double_t SumSNRD=0, SumSNRR=0,SumSNR=0;
  int chanD[2]={0,0};
  int chanR[2]={0,0};
  int chanDsame=0,chanRsame=0,chanRsameB=0;
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(p[2*TotalAntennasRx +iRx]!=0){
      chanD[0]++;
    }
    if(p[3*TotalAntennasRx +iRx]!=0){
      chanR[0]++;
    }
    if(IgnoreCh[0][iRx]!=0){
      chanD[1]++;
    }
    if(IgnoreCh[1][iRx]!=0){
      chanR[1]++;
    }
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
      SumSNRD+=p[4*TotalAntennasRx +iRx];
      chanDsame++;
    }
    if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      SumSNRR+=p[5*TotalAntennasRx +iRx];
      chanRsame++;
    }
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      chanRsameB++;
    }

  }
 
  double deltaT=0;
  double output1=0;  
  double output2=0;

  //if((p[6*TotalAntennasRx+13]==1 && chanDsame>=chanD[0]) || ( p[6*TotalAntennasRx+13]==0 && chanRsameB>=chanD[0])){
    if(chanDsame>=chanD[0]){

    //if(p[6*TotalAntennasRx+13]==1){
      double loop1_D=0,loop1_R=0;   
      double loop2_D=0,loop2_R=0;
      double sum1_D=0,sum1_R=0;  

      for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
	//cout<<"ray times are "<<iRx<<" "<<p[0+iRx]<<" "<<p[4*TotalAntennasRx +iRx]<<" "<<p[1*TotalAntennasRx+iRx]<<" "<<p[5*TotalAntennasRx +iRx]<<endl;
	//cout<<"sim ray times are "<<iRx<<" "<<timeRay[0][iRx]<<" "<<timeRay[1][iRx]<<endl;
	loop2_D=0,loop2_R=0;
	double sum_D=0,sum_R=0;
	double sum_W_D=0,sum_W_R=0;

	for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 
	  if(p[2*TotalAntennasRx +iRx]!=0 && p[2*TotalAntennasRx +iRxB]!=0 && IgnoreCh[0][iRx]!=0 && IgnoreCh[0][iRxB]!=0){     

	    if(sum_D==0){
	      loop1_D+=p[4*TotalAntennasRx +iRx];
	    }
	    deltaT=(timeRay[0][iRx] - timeRay[0][iRxB] - (p[0+iRx] - p[0+iRxB]));
	      
	    //cout<<"delta D are  "<<iRx<<" "<<iRxB<<" "<<deltaT<<" "<<timeRay[0][iRx]<<" "<<timeRay[0][iRxB]<<" "<<p[0+iRx]<<" "<<p[0+iRxB]<<endl;
	    
	    sum_D+=pow(deltaT,2)*p[4*TotalAntennasRx +iRxB];
	    sum_W_D+=p[4*TotalAntennasRx +iRxB];
	    loop2_D++;

	  }
	}

	if(loop2_D>0){
	  sum_D=sum_D/sum_W_D;
	  if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
	    sum1_D+=sum_D*p[4*TotalAntennasRx +iRx];
	  }
	}

	loop2_R=0;
	sum_R=0;
	sum_W_R=0;

	for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 
	  if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0  && p[3*TotalAntennasRx +iRxB]!=0 && IgnoreCh[1][iRxB]!=0){
	    if(sum_R==0){
	      loop1_D+=p[5*TotalAntennasRx +iRx];
	      sum_R=1;
	    }
	    deltaT=(timeRay[1][iRx] - timeRay[1][iRxB] - (p[1*TotalAntennasRx+iRx] - p[1*TotalAntennasRx+iRxB]));
	    //cout<<"delta R are  "<<iRx<<" "<<iRxB<<" "<<deltaT<<" "<<timeRay[1][iRx]<<" "<<timeRay[1][iRxB]<<" "<<p[1*TotalAntennasRx+iRx]<<" "<<p[1*TotalAntennasRx+iRxB]<<endl;
	    sum_D+=pow(deltaT,2)*p[5*TotalAntennasRx +iRxB];
	    sum_W_D+=p[5*TotalAntennasRx +iRxB];
	    loop2_R++;
	  }
	}

	if(loop2_R>0){
	  sum_R=sum_D/sum_W_D;
	  if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
	    sum1_D+=sum_R*p[5*TotalAntennasRx +iRx];
	  }
	}

      }
      
      output1+=sum1_D/loop1_D;
    // }

    // if(p[6*TotalAntennasRx+13]==0){
    //   double loop1=0;
    //   double loop2=0;
    
    //   double sum1=0;  
    //   for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
    // 	loop2=0;
    // 	double sum=0;
    // 	double sum_W=0;
    // 	for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 

    // 	  if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0 && p[2*TotalAntennasRx +iRxB]!=0 && IgnoreCh[1][iRxB]!=0){

    // 	    if(sum==0){
    // 	      loop1+=p[4*TotalAntennasRx +iRx];
    // 	    }

    // 	    deltaT=(timeRay[1][iRx] - timeRay[1][iRxB]  - (p[0+iRx] - p[0+iRxB]));
    // 	    sum+=pow(deltaT,2)*p[4*TotalAntennasRx +iRxB];
    // 	    sum_W+=p[4*TotalAntennasRx +iRxB];	  
    // 	    loop2++;

    // 	  }
    // 	}
      
    // 	if(loop2>0){
    // 	  sum=sum/sum_W;
    // 	  if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
    // 	    sum1+=sum*p[4*TotalAntennasRx +iRx];
    // 	  }
    // 	}
    //   }
      
    //   output2+=sum1/loop1;

    // }
	  
    
      //if(p[6*TotalAntennasRx+13]==1){
      output=output1;
      //}
 
      //if(p[6*TotalAntennasRx+13]==0){
      //output=output2;
      //}  

    }else{
    output=1e9+r;
  }

  if(std::isnan(output)){
    output=1e9+r;
  }
  
  if((gsl_vector_get(v, 0)>179.9 || gsl_vector_get(v, 0)<0.1)){
    output=1e9+output;
  }

  if(XYZ[2]+AvgAntennaCoordRx[2]>0 && output<1e9){
    output=1e9+output;
  }

  return output;
}

double Interferometer::Minimizer_fCnz(const gsl_vector *v, void *params){

  double *p = (double *)params; 
  
  double ExpectedX=p[6*TotalAntennasRx];
  double ExpectedY=p[6*TotalAntennasRx+1];
  double ExpectedZ=p[6*TotalAntennasRx+2];

  double ExpectedTh=p[6*TotalAntennasRx+3];
  double ExpectedPh=p[6*TotalAntennasRx+4];
  double ExpectedR=p[6*TotalAntennasRx+5];

  double theta, phi, r;
  theta = gsl_vector_get(v, 0)*(Interferometer::pi/180.0);
  phi = gsl_vector_get(v, 1)*(Interferometer::pi/180.0);
  r = p[6*TotalAntennasRx+10];

  Double_t ThPhR[3]={theta,phi,r};
  Double_t XYZ[3]={0,0,0};
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
  
  double output=0;
  
  double AntennaCoordTx[3]={0,0,0};
  for(int ixyz=0;ixyz<3;ixyz++){
    AntennaCoordTx[ixyz]=XYZ[ixyz];
  }
    
  vector<double> timeRay[2];
  vector<int> IgnoreCh[2];

  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }

  if(XYZ[2]+AvgAntennaCoordRx[2]<0 || p[6*TotalAntennasRx+12]==0){
    //cout<<"working for ice "<<gsl_vector_get(v, 0)<<" "<<gsl_vector_get(v, 1)<<" "<<p[6*TotalAntennasRx+10]<<" "<<XYZ[2]+AvgAntennaCoordRx[2]<<" "<<p[6*TotalAntennasRx+12]<<endl;;
    Interferometer::GenerateChHitTimeAndCheckHits_Cnz(AntennaCoordTx,timeRay,IgnoreCh);  
  }else{
    //cout<<"working for air "<<gsl_vector_get(v, 0)<<" "<<gsl_vector_get(v, 1)<<" "<<p[6*TotalAntennasRx+10]<<" "<<XYZ[2]+AvgAntennaCoordRx[2]<<" "<<p[6*TotalAntennasRx+12]<<endl;;
    Interferometer::GenerateChHitTimeAndCheckHits_Air(AntennaCoordTx,timeRay,IgnoreCh);  
  }
  
  Double_t SumSNRD=0, SumSNRR=0,SumSNR=0;
  int chanD[2]={0,0};
  int chanR[2]={0,0};
  int chanDsame=0,chanRsame=0,chanRsameB=0;
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(p[2*TotalAntennasRx +iRx]!=0){
      chanD[0]++;
    }
    if(p[3*TotalAntennasRx +iRx]!=0){
      chanR[0]++;
    }
    if(IgnoreCh[0][iRx]!=0){
      chanD[1]++;
    }
    if(IgnoreCh[1][iRx]!=0){
      chanR[1]++;
    }
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
      SumSNRD+=p[4*TotalAntennasRx +iRx];
      chanDsame++;
    }
    if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      SumSNRR+=p[5*TotalAntennasRx +iRx];
      chanRsame++;
    }
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      chanRsameB++;
    }

  }
 
  double deltaT=0;
  double output1=0;  
  double output2=0;

  // if((p[6*TotalAntennasRx+13]==1 && chanDsame>=chanD[0]) || ( p[6*TotalAntennasRx+13]==0 && chanRsameB>=chanD[0] && chanDsame>=chanD[0])){
  if(chanDsame>=chanD[0]){
    //if(p[6*TotalAntennasRx+13]==1){

    double loop1_D=0,loop1_R=0;   
    double loop2_D=0,loop2_R=0;
    double sum1_D=0,sum1_R=0;  

    for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
     
      loop2_D=0,loop2_R=0;
      double sum_D=0,sum_R=0;
      double sum_W_D=0,sum_W_R=0;

      for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 
	if(p[2*TotalAntennasRx +iRx]!=0 && p[2*TotalAntennasRx +iRxB]!=0 && IgnoreCh[0][iRx]!=0 && IgnoreCh[0][iRxB]!=0){     

	  if(sum_D==0){
	    loop1_D+=p[4*TotalAntennasRx +iRx];
	  }
	    
	  deltaT=(timeRay[0][iRx] - timeRay[0][iRxB] - (p[0+iRx] - p[0+iRxB]));
	  sum_D+=pow(deltaT,2)*p[4*TotalAntennasRx +iRxB];
	  sum_W_D+=p[4*TotalAntennasRx +iRxB];
	  loop2_D++;

	}
      }

      if(loop2_D>0){
	sum_D=sum_D/sum_W_D;
	if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
	  sum1_D+=sum_D*p[4*TotalAntennasRx +iRx];
	}
      }

      // loop2_R=0;
      // sum_R=0;
      // sum_W_R=0;

      // for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 
      //   if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0  && p[3*TotalAntennasRx +iRxB]!=0 && IgnoreCh[1][iRxB]!=0){
      //     if(sum_R==0){
      // 	loop1_D+=p[5*TotalAntennasRx +iRx];
      // 	sum_R=1;
      //     }

      //     deltaT=(timeRay[1][iRx] - timeRay[1][iRxB] - (p[1*TotalAntennasRx+iRx] - p[1*TotalAntennasRx+iRxB]));
      //     sum_D+=pow(deltaT,2)*p[5*TotalAntennasRx +iRxB];
      //     sum_W_D+=p[5*TotalAntennasRx +iRxB];
      //     loop2_R++;
      //   }
      // }

      // if(loop2_R>0){
      //   sum_R=sum_D/sum_W_D;
      //   if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      //     sum1_D+=sum_R*p[5*TotalAntennasRx +iRx];
      //   }
      // }

    }
	      
    output1+=sum1_D/loop1_D;
    // if(loop1_R>0){
    //   output1+=sum1_R/loop1_R;
    // }
    //}

    // if(p[6*TotalAntennasRx+13]==0){
    // 	double loop1=0;
    // 	double loop2=0;
    
    // 	double sum1=0;  
    // 	for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
    // 	  loop2=0;
    // 	  double sum=0;
    // 	  double sum_W=0;
    // 	  for(int iRxB=0;iRxB<TotalAntennasRx;iRxB++){ 

    // 	    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0 && p[2*TotalAntennasRx +iRxB]!=0 && IgnoreCh[1][iRxB]!=0){

    // 	      if(sum==0){
    // 		loop1+=p[4*TotalAntennasRx +iRx];
    // 	      }

    // 	      deltaT=(timeRay[1][iRx] - timeRay[1][iRxB]  - (p[0+iRx] - p[0+iRxB]));
    // 	      sum+=pow(deltaT,2)*p[4*TotalAntennasRx +iRxB];
    // 	      sum_W+=p[4*TotalAntennasRx +iRxB];	  
    // 	      loop2++;

    // 	    }
    // 	  }
      
    // 	  if(loop2>0){
    // 	    sum=sum/sum_W;
    // 	    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
    // 	      sum1+=sum*p[4*TotalAntennasRx +iRx];
    // 	    }
    // 	  }
    // 	}
      
    // 	output2+=sum1/loop1;

    // }

    //if(p[6*TotalAntennasRx+13]==1){
    output=output1;
    //}
 
    // if(p[6*TotalAntennasRx+13]==0){
    //   output=output2;
    // }
    

  }else{
    output=1e9+r;
  }

  if(std::isnan(output)){
    output=1e9+r;
  }
  
  if((gsl_vector_get(v, 0)>179.9 || gsl_vector_get(v, 0)<0.1)){
    output=1e9+output;
  }

  if(XYZ[2]+AvgAntennaCoordRx[2]>0 && output<1e9){
    output=1e9+output;
  }

  return output;
}

void Interferometer::MinimizerThPhR1D(double m, double a, double b, double &FinalMinValue, double &FnMinValue, int &Iterations, void *parameters){

  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;

  F.function = &Interferometer::Minimizer_f1D;
  F.params = parameters;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  // printf ("using %s method\n",
  //         gsl_min_fminimizer_name (s));

  // printf ("%5s [%9s, %9s] %9s %9s\n",
  //         "iter", "lower", "upper", "min", "err(est)");

  // printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
  //         iter, a, b,
  //         m, b - a);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status
        = gsl_min_test_interval (a, b, 0.001, 0.0);

      // if (status == GSL_SUCCESS)
      //   printf ("Converged:\n");

      // printf ("%5d [%.7f, %.7f] "
      //         "%.7f %.7f\n",
      //         iter, a, b,
      //         m, b - a);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  Iterations=iter;
  FinalMinValue=m;
  FnMinValue=(*((F).function))(m,(F).params);

  gsl_min_fminimizer_free (s);
}

double Interferometer::MinimizerThPh(double x, void * params){
 
  double *p = (double *)params;

  double ThPhR[3]={p[6*TotalAntennasRx+3]*(Interferometer::pi/180.),p[6*TotalAntennasRx+4]*(Interferometer::pi/180.),x};
  double TxCor[3]={0,0,0};
  Interferometer::ThPhRtoXYZ(ThPhR,TxCor);
  
  vector<double> timeRay[2];
  vector<int> IgnoreCh[2];

  timeRay[0].resize(TotalAntennasRx);
  timeRay[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
  }
  Interferometer::GenerateChHitTimeAndCheckHits(TxCor,timeRay,IgnoreCh);

  double FinalMinValue=0;
  p[6*TotalAntennasRx+10]=x;

  const gsl_multimin_fminimizer_type *MinimisationType = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *MinimizerWorkSpace = NULL;
  gsl_vector *XYZStepSizeVector, *XYZVec;
  gsl_multimin_function MultiDimMinFunc;

  size_t Iter = 0;
  int Status;
  double Size;
    
  /* Starting point */
  XYZVec = gsl_vector_alloc (2);
  gsl_vector_set (XYZVec, 0, p[6*TotalAntennasRx+3]);
  gsl_vector_set (XYZVec, 1, p[6*TotalAntennasRx+4]);
  
  double ExpectedDisplacement=p[6*TotalAntennasRx+5];
  double ThetaPhiStepSize=atan(Interferometer::ExpectedUncertainty/ExpectedDisplacement)*(180./Interferometer::pi);
  
  /* Set initial step sizes */
  XYZStepSizeVector = gsl_vector_alloc (2);
  gsl_vector_set (XYZStepSizeVector, 0, ThetaPhiStepSize*4);
  gsl_vector_set (XYZStepSizeVector, 1, ThetaPhiStepSize*4);
  
  /* Initialize method and iterate */
  MultiDimMinFunc.n = 2;
  MultiDimMinFunc.f = Minimizer_f;
  MultiDimMinFunc.params = p;

  MinimizerWorkSpace = gsl_multimin_fminimizer_alloc (MinimisationType, 2);
  gsl_multimin_fminimizer_set (MinimizerWorkSpace, &MultiDimMinFunc, XYZVec, XYZStepSizeVector);


  double FinalTxCor[2]; 
  int Iterations=0;
  int Max_Iter=p[6*TotalAntennasRx+14];
  do
    {
      Iter++;
      Iterations++;
      Status = gsl_multimin_fminimizer_iterate(MinimizerWorkSpace);
     
      if (Status){
        break;
      }

      Size = gsl_multimin_fminimizer_size (MinimizerWorkSpace);
      Status = gsl_multimin_test_size (Size, 1e-3);
      
      // if (Status == GSL_SUCCESS)
      //   {
      //     printf ("converged to minimum at\n");
      //   }
      // printf ("%5d %7.3f %7.3f %7.3f f() = %7.3f size = %.3f\n",
      //         Iter,
      //         gsl_vector_get (MinimizerWorkSpace->x, 0),
      //         gsl_vector_get (MinimizerWorkSpace->x, 1),
      // 	      x,
      //         MinimizerWorkSpace->fval,Size);
     
      // cout<<Iter<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 0)<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 1) <<" "<<MinimizerWorkSpace->fval<<" "<<Size<<" "<<(*((MultiDimMinFunc).f))(XYZVec,(MultiDimMinFunc).params)<<endl;
      
      FinalMinValue=MinimizerWorkSpace->fval;
      FinalTxCor[0]=gsl_vector_get (MinimizerWorkSpace->x, 0);
      FinalTxCor[1]=gsl_vector_get (MinimizerWorkSpace->x, 1);
    }
  while (Status == GSL_CONTINUE && Iter < Max_Iter);

  //printf ("error: %s\n", gsl_strerror (Status));
  p[6*TotalAntennasRx+11]=FinalMinValue;
  
  p[6*TotalAntennasRx+6]=FinalTxCor[0];
  p[6*TotalAntennasRx+7]=FinalTxCor[1];
  p[6*TotalAntennasRx+8]=x;

  gsl_vector_free(XYZVec);
  gsl_vector_free(XYZStepSizeVector);
  gsl_multimin_fminimizer_free (MinimizerWorkSpace);
  
  return FinalMinValue;
}

void Interferometer::MinimizerThPhR(double InitialTxCor_XYZ[3], double InitialTxCor_ThPhR[3], double FinalTxCor[3], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], double &FinalMinValue, int &Iterations, double MinimizerRadialWidth, int IsItBelowStation, int max_iter){
 
  int wloop=0;
  double m = InitialTxCor_ThPhR[2]; 
  double a = InitialTxCor_ThPhR[2]-MinimizerRadialWidth, b = InitialTxCor_ThPhR[2]+MinimizerRadialWidth;
  m=a;

  bool HitBoundary=true;

  double ParameterArray[6*TotalAntennasRx+15];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }

  ParameterArray[iEnt]=InitialTxCor_XYZ[0];
  ParameterArray[iEnt+1]=InitialTxCor_XYZ[1];
  ParameterArray[iEnt+2]=InitialTxCor_XYZ[2];
  
  ParameterArray[iEnt+3]=InitialTxCor_ThPhR[0];
  ParameterArray[iEnt+4]=InitialTxCor_ThPhR[1];
  ParameterArray[iEnt+5]=InitialTxCor_ThPhR[2];

  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=0;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;

  if(InitialTxCor_ThPhR[0]<45){
    ParameterArray[iEnt+12]=1;
  }
  
  ParameterArray[iEnt+13]=IsItBelowStation;
  ParameterArray[iEnt+14]=max_iter;

  m = InitialTxCor_ThPhR[2]; 
  a = InitialTxCor_ThPhR[2]-MinimizerRadialWidth, b = InitialTxCor_ThPhR[2]+MinimizerRadialWidth;
    
  int status;
  int iter = 0, Max_Iter = 50;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s; 
 
  if(a<40){
    a=40;
  }
 
  if(fabs(m-40)<1){
    m=(a+b)/2;
  }
  double a_ini=a;
  double b_ini=b;
  //cout<<"Minimizer input values are  a "<<a<<", m "<<m<<", b "<<b<<endl;
  //cout<<InitialTxCor_ThPhR[0]<<" "<<InitialTxCor_ThPhR[1]<<" "<<InitialTxCor_ThPhR[2]<<" "<<endl;
  gsl_function F;
  
  F.function = &MinimizerThPh;   
  F.params = ParameterArray;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);
 
  // printf ("using %s method\n",
  //         gsl_min_fminimizer_name (s));

  // printf ("%5s [%9s, %9s] %9s %10s %9s\n",
  //         "iter", "lower", "upper", "min",
  //         "err", "err(est)");

  // printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
  //         iter, a, b,
  //         m, b - a);  
  
  Iterations=0;
  do
    {
      
      iter++;
      Iterations++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      a = gsl_min_fminimizer_x_lower (s);
      b = gsl_min_fminimizer_x_upper (s);

      status= gsl_min_test_interval (a, b, 1e-2,0);

      // if (status == GSL_SUCCESS)
      //   printf ("Converged:\n");
      
      // printf ("%5d [%.7f, %.7f] "
      //         "%.7f %.7f\n",
      //         iter, a, b,
      //         m, b - a);

      // if(fabs(m-a_ini)<10 || fabs(m-b_ini)<10){
      //   cout<<"Minimizer loop stopped as R value came to close to one of the boundary values. Iterations: "<<Iterations<<endl;
      //   status = GSL_SUCCESS;
      // }

    }
  while (status == GSL_CONTINUE && iter < Max_Iter);
  
  gsl_min_fminimizer_free (s);
  FinalTxCor[0]=ParameterArray[iEnt+6];
  FinalTxCor[1]=ParameterArray[iEnt+7];
  FinalTxCor[2]=m;
  FinalMinValue=ParameterArray[6*TotalAntennasRx+11]; 
  b=b_ini;
  a=a_ini;
  HitBoundary=false;
  
  if(fabs(m-b_ini)<10 || fabs(m-a_ini)<10){
    HitBoundary=true;     
    cout<<"Entering secondary minimization with reconstructed values: Th "<<FinalTxCor[0]<<", Ph "<<FinalTxCor[1]<<", R "<<FinalTxCor[2]<<", FnMin "<<FinalMinValue<<endl;
  
    double StartR=FinalTxCor[2]-MinimizerRadialWidth;   
    double StopR=5000;
    double StepSizeR=50;

    if(fabs(m-b_ini)<10){
      StartR=FinalTxCor[2]-MinimizerRadialWidth*2;   
      //StopR=FinalTxCor[2]+MinimizerRadialWidth/2;
      StopR=5000;
      // if(b_ini>=4000-10){
      // 	StopR=500;
      // }

    }

    if(fabs(m-a_ini)<10){
      //StartR=FinalTxCor[2]-MinimizerRadialWidth/2;   
      StartR=40;   
      StopR=FinalTxCor[2]+MinimizerRadialWidth*2;
    }

    if(StopR<1000){
      StopR=1000;
    }
      
    if(StartR<40){
      StartR=40;
    }

    double testTh=FinalTxCor[0];
    double testPh=FinalTxCor[1];
    double testR=0;

    int Iterations_1D=0;
    double FnMinValue_1D;

    ParameterArray[iEnt+3]=FinalTxCor[0];
    ParameterArray[iEnt+4]=FinalTxCor[1];
    ParameterArray[iEnt+5]=FinalTxCor[2];

    Interferometer::MinimizerThPhR1D((StartR+StopR)/2, StartR, StopR, FinalTxCor[2], FnMinValue_1D, Iterations_1D, ParameterArray);

    Iterations=Iterations_1D;
    FinalMinValue=FnMinValue_1D; 

  }

}

double Interferometer::GetChiSquaredThPhR(double UserCor[3], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], int IsItBelowStation,int CnzOrEnz, int max_iter){

  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  double ParameterArray[6*TotalAntennasRx+15];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=0;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  ParameterArray[iEnt+13]=IsItBelowStation;

  ParameterArray[iEnt+14]=max_iter;

  double Tht=UserCor[0];
  double Pht=UserCor[1];
  double Rt=UserCor[2];
    
  if(Tht<45){
    ParameterArray[iEnt+12]=1;
  }

  Double_t ThPhR[3]={Tht*(Interferometer::pi/180),Pht*(Interferometer::pi/180),Rt};
  Double_t XYZ[3]={0,0,0}; 
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);

  ParameterArray[iEnt]=XYZ[0];
  ParameterArray[iEnt+1]=XYZ[1];
  ParameterArray[iEnt+2]=XYZ[2];
  ParameterArray[iEnt+3]=Tht;
  ParameterArray[iEnt+4]=Pht;
  ParameterArray[iEnt+5]=Rt;      

  gsl_vector_set (ThPhRVec, 0, Tht);
  gsl_vector_set (ThPhRVec, 1, Pht);
  ParameterArray[iEnt+10]=Rt;

  double min=0;
  if(CnzOrEnz==0){
    min=Interferometer::Minimizer_f(ThPhRVec, ParameterArray);
  }
  if(CnzOrEnz==1){
    min=Interferometer::Minimizer_fCnz(ThPhRVec, ParameterArray);
  }

  gsl_vector_free(ThPhRVec);

  return min;
  
}

void Interferometer::GetRecieveAngle(vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], int max_iter, double StartCor[3]){

  double TestCor[3]={0,0,0}; 
  double min=1e9;
  int minbin=0;

  vector <int> IgnoreChB[2];
  vector <double> ChSNRB[2];
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    IgnoreChB[0].push_back(IgnoreCh[0][iRx]);
    ChSNRB[0].push_back(ChSNR[0][iRx]);
    IgnoreChB[1].push_back(0);
    ChSNRB[1].push_back(1);
  }
  
  for(int ith=1;ith<90;ith++){
    for(int iph=0;iph<180;iph++){
      TestCor[0]=2*ith;   
      TestCor[1]=2*iph-180;  
      TestCor[2]=5;
  
      double X2a=Interferometer::GetChiSquaredThPhR(TestCor, ChHitTime, IgnoreChB, ChSNRB, 1,1, max_iter);
   
      if(std::isnan(X2a)==false && X2a<1e9){   
	if(X2a<min){
	  min=X2a;
	  minbin=ith;

	  StartCor[0]=TestCor[0];
	  StartCor[1]=TestCor[1];
	  StartCor[2]=TestCor[2];
	}
      }
    }
  }

  cout<<"Approximate Arrival Direction is: Theta= "<<StartCor[0]<<" , Phi= "<<StartCor[1]<<endl;//" "<<StartCor[2]<<" "<<min<<endl;
}

void Interferometer::SearchApproxiMin(double GuessResultCor[3][4],double ParameterArray[6*TotalAntennasRx+15], bool CheckBelowSurface){
  
  double StartCor[3];

  int TotEnt=6*TotalAntennasRx;

  vector<int> IgnoreCh[2];
  vector<double> ChHitTime[2];
  vector <double> ChSNR[2];

  ChHitTime[0].resize(TotalAntennasRx);
  ChHitTime[1].resize(TotalAntennasRx);
  IgnoreCh[0].resize(TotalAntennasRx);
  IgnoreCh[1].resize(TotalAntennasRx);
  ChSNR[0].resize(TotalAntennasRx);
  ChSNR[1].resize(TotalAntennasRx);
  
  int iEnt=0;

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChHitTime[iray][iRx]=ParameterArray[iEnt];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      IgnoreCh[iray][iRx]=ParameterArray[iEnt];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ChSNR[iray][iRx]=ParameterArray[iEnt];
      iEnt++;
    }
  }

  double min=1e9;
  int minbin=0;
  double PossibleMin[2]={0,0};
  double PossibleMinCor[2][3]={{0,0,0},{0,0,0}};

  min=1e9;
  minbin=0;

  double TestCor1[3]={95,-180,250};
   
  for(int idir=0;idir<70;idir++){
     TestCor1[1]=-180+5*idir;
     double X2a=Interferometer::GetChiSquaredThPhR(TestCor1 , ChHitTime, IgnoreCh, ChSNR, ParameterArray[TotEnt+13],1,ParameterArray[TotEnt+14]);

     if(std::isnan(X2a)==false && X2a<1e9){      
       if(X2a<min){
	 min=X2a;
	 minbin=idir;

	 StartCor[0]=TestCor1[0];
	 StartCor[1]=TestCor1[1];
	 StartCor[2]=TestCor1[2];
       }
     }
   }

  double TestCor[3]={0,StartCor[1],0};   
  min=1e9;
  minbin=0;

  int StartThBin=1;
  if(CheckBelowSurface==true){
    //cout<<"checking below surface"<<endl;  
    StartThBin=18;
  }

  for(int idir=1;idir<72;idir++){
    for(int idis=4;idis<50;idis++){
    //for(int idir=32;idir<72;idir++){

      TestCor[0]=2.5*idir;   
      TestCor[2]=50*idis;  
      
      double X2a=Interferometer::GetChiSquaredThPhR(TestCor , ChHitTime, IgnoreCh, ChSNR, ParameterArray[TotEnt+13],0, ParameterArray[TotEnt+14]);
    
      if(std::isnan(X2a)==false && X2a<1e9){   
	//cout<<TestCor[0]<<" "<<TestCor[2]<<" "<<X2a<<endl;
	if(X2a<min){
	  min=X2a;
	  minbin=idir;

	  StartCor[0]=TestCor[0];
	  StartCor[1]=TestCor[1];
	  StartCor[2]=TestCor[2];
	}
      }
    }
  }

  cout<<" minimum coordinates are "<<StartCor[0]<<" "<<StartCor[1]<<" "<<StartCor[2]<<" "<<min<<endl;
 
  Double_t NumBinsTh=20,NumBinsPh=10,NumBinsR=5;
  Double_t StartTh=StartCor[0]-15,StartPh=StartCor[1]-10,StartR=200;
  Double_t StopTh=StartCor[0]+15,StopPh=StartCor[1]+10,StopR=2500;
 
  if(StartTh<=2){
    StartTh=2;
  } 
  if(StopTh>=178){
    StopTh=178;
  }
  if(StartR<40){
    StartR=40;
  }

  // cout<<"start and stop thetas are "<<StartTh<<" "<<StopTh<<endl;
  // cout<<"start and stop phis are "<<StartPh<<" "<<StopPh<<endl;
  // cout<<"start and stop Rs are "<<StartR<<" "<<StopR<<endl;

  Double_t StepSizeTh=(StopTh-StartTh)/NumBinsTh,StepSizePh=(StopPh-StartPh)/NumBinsPh;
  Double_t StepSizeR=(StopR-StartR)/NumBinsR;
  
  vector <double> RecoPar[4];
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);

  for(double i=0; i<=NumBinsTh;i++){
    double Tht=i*StepSizeTh + StartTh;
    for(double j=0; j<=NumBinsPh;j++){	
      double Pht=j*StepSizePh + StartPh;
      for(double k=0; k<=NumBinsR;k++){
  	double Rt=k*StepSizeR + StartR;
	
	Double_t ThPhRa[3]={Tht,Pht,Rt};
	double min=Interferometer::GetChiSquaredThPhR(ThPhRa, ChHitTime, IgnoreCh, ChSNR, ParameterArray[TotEnt+13],0, ParameterArray[TotEnt+14]);

  	if(std::isnan(min)==false && min!=1e9 && min<1e9 && min!=0){
  	  RecoPar[0].push_back(Tht);
  	  RecoPar[1].push_back(Pht);
  	  RecoPar[2].push_back(Rt);
  	  RecoPar[3].push_back(min);
  	}
      }
    }
  }
  gsl_vector_free(ThPhRVec);
  
  int Nmin=0;
  int Nminsize=3;
  if(RecoPar[3].size()<3){
    Nminsize=RecoPar[3].size();
  }
  while(Nmin<Nminsize){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=RecoPar[2][FinalMinValueBin];
    GuessResultCor[Nmin][3]=RecoPar[3][FinalMinValueBin];
    //cout<<"Guess Minimas: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }

  cout<<"First guess at Theta and Phi: Tht "<<GuessResultCor[0][0]<<", Pht "<<GuessResultCor[0][1]<<" R "<<GuessResultCor[0][2]<<" with min "<<GuessResultCor[0][3]<<endl;
  
}

void Interferometer::GetApproximateMinUserCor(vector <double>(&UserCor)[3],double GuessResultCor[3][4], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], int IsItBelowStation, int max_iter){
  
  vector <double> RecoPar[4];
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  double ParameterArray[6*TotalAntennasRx+15];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=0;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  ParameterArray[iEnt+13]=IsItBelowStation;

  ParameterArray[iEnt+14]=max_iter;

  for(int i=0;i<UserCor[0].size();i++){
    double Tht=UserCor[0][i];
    double Pht=UserCor[1][i];
    double Rt=UserCor[2][i];
    
    if(Tht<45){
      ParameterArray[iEnt+12]=1;
    }else{
      ParameterArray[iEnt+12]=0;
    }
    
    Double_t ThPhR[3]={Tht*(Interferometer::pi/180),Pht*(Interferometer::pi/180),Rt};
    Double_t XYZ[3]={0,0,0}; 
    Interferometer::ThPhRtoXYZ(ThPhR,XYZ);

    ParameterArray[iEnt]=XYZ[0];
    ParameterArray[iEnt+1]=XYZ[1];
    ParameterArray[iEnt+2]=XYZ[2];
    ParameterArray[iEnt+3]=Tht;
    ParameterArray[iEnt+4]=Pht;
    ParameterArray[iEnt+5]=Rt;      

    gsl_vector_set (ThPhRVec, 0, Tht);
    gsl_vector_set (ThPhRVec, 1, Pht);
    ParameterArray[iEnt+10]=Rt;
   
    double min=0;
    min=Interferometer::Minimizer_f(ThPhRVec, ParameterArray);
   
    if(std::isnan(min)==false && min!=1e9 && min<1e9 && Tht>0.001){
      
      RecoPar[0].push_back(Tht);
      RecoPar[1].push_back(Pht);
      RecoPar[2].push_back(Rt);
      RecoPar[3].push_back(min);
    }
  }
  gsl_vector_free(ThPhRVec);
  
  int Nmin=0;
  int Nminsize=2;
  if(RecoPar[3].size()<2){
    Nminsize=RecoPar[3].size();
  }
  while(Nmin<Nminsize){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=RecoPar[2][FinalMinValueBin];
    GuessResultCor[Nmin][3]=RecoPar[3][FinalMinValueBin];
    cout<<"Guess Minimas: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }
}

void Interferometer::GetApproximateMinThPhR(double GuessResultCor[3][4], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], int IsItBelowStation, int max_iter){
  
  double ParameterArray[6*TotalAntennasRx+15];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=0;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;

  ParameterArray[iEnt+14]=max_iter;

  bool CheckBelowSurface=false;
 
  ParameterArray[iEnt+13]=IsItBelowStation;

  if(IsItBelowStation==0){
    ParameterArray[iEnt+13]=0;

    bool AreAllRIgnored=true;
    bool IsdtDRtimeInWindow=false;
    
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(IgnoreCh[1][iRx]==1){
	double dtDR=fabs(ChHitTime[1][iRx]-ChHitTime[0][iRx]);
	if(dtDR>500){
	  AreAllRIgnored=true;
	  IsdtDRtimeInWindow=false;
	}else{
	  AreAllRIgnored=false;
	  IsdtDRtimeInWindow=true;
	}
      }
    }
    if(AreAllRIgnored==true && IsdtDRtimeInWindow==false){
      CheckBelowSurface=false;
      cout<<"The vertex could be above or below station. No R rays were detected."<<endl;
      ParameterArray[iEnt+12]=1;
      Interferometer::SearchApproxiMin(GuessResultCor,ParameterArray,CheckBelowSurface);
    }

    if(AreAllRIgnored==false && IsdtDRtimeInWindow==true){
      CheckBelowSurface=true;
      cout<<"The vertex is probably below station. Some R rays were detected."<<endl;
      ParameterArray[iEnt+12]=0;
      Interferometer::SearchApproxiMin(GuessResultCor,ParameterArray,CheckBelowSurface);
    }
    
  }

  if(IsItBelowStation==1){
    cout<<"The vertex is probably below station"<<endl;
    CheckBelowSurface=true;
    ParameterArray[iEnt+12]=0;
    ParameterArray[iEnt+13]=1;
    Interferometer::SearchApproxiMin(GuessResultCor,ParameterArray,CheckBelowSurface);
  }

}

void Interferometer::GetApproximateDistance(double GuessResultCor[3][4], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], int IsItBelowStation, int max_iter){
  
  double ParameterArray[6*TotalAntennasRx+15];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=0;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  ParameterArray[iEnt+13]=IsItBelowStation;

  if(GuessResultCor[0][0]<45){
    ParameterArray[iEnt+12]=1;
  }

  ParameterArray[iEnt+14]=max_iter;

  double min=0;
  double GuessFinalMin[3]={0,0,0};
  vector <double> RecoPar[4];

  double StartR=3000;
  double StopR=200;

  if(StopR<40){
    StopR=40;
  }  

  double StepSizeR=50;//(StopR-StartR)/20;
  double FinalTxCor[3];
  bool checkmincurve=false;
  int iloop=0;
  
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  int countnan=0;
  
  double testTh=GuessResultCor[0][0];
  double testPh=GuessResultCor[0][1];
  double testR=StopR+1;
  
  while(checkmincurve==false && testR>StopR){
    
    testTh=GuessResultCor[0][0];
    testPh=GuessResultCor[0][1];
    testR=StartR - (iloop)*StepSizeR ;
      
    Double_t ThPhR[3]={testTh*(Interferometer::pi/180),testPh*(Interferometer::pi/180),testR};
    Double_t XYZ[3]={0,0,0}; 
    Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
    ParameterArray[iEnt]=XYZ[0];
    ParameterArray[iEnt+1]=XYZ[1];
    ParameterArray[iEnt+2]=XYZ[2];
    ParameterArray[iEnt+3]=testTh;
    ParameterArray[iEnt+4]=testPh;
    ParameterArray[iEnt+5]=testR;

    gsl_vector_set (ThPhRVec, 0, testTh);
    gsl_vector_set (ThPhRVec, 1, testPh);
    ParameterArray[iEnt+10]=testR;
      
    min=Interferometer::Minimizer_f(ThPhRVec, ParameterArray);

    FinalTxCor[0]=testTh;
    FinalTxCor[1]=testPh;
    FinalTxCor[2]=testR;
    
    if(std::isnan(min)==false && min!=1e9 && min<1e9 && FinalTxCor[0]>0.001 && min!=0){
      RecoPar[0].push_back(FinalTxCor[0]);
      RecoPar[1].push_back(FinalTxCor[1]);
      RecoPar[2].push_back(FinalTxCor[2]);
      RecoPar[3].push_back(min);
    }
    
    int size=RecoPar[3].size();
    if(size>0 && iloop>=5){
      //if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3]) ||  (fabs(RecoPar[3][size-1]-RecoPar[3][size-2])/RecoPar[3][size-1]<0.01 && fabs(RecoPar[3][size-3]-RecoPar[3][size-3])/RecoPar[3][size-2]<0.01) ){
      // if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3] && iloop>=5) ){
      //   cout<<"found the minima!"<<endl;
      //   checkmincurve=true;
      // }
    }
    // if((std::isnan(min)==true || min==1e9 || min>1e9 || FinalTxCor[0]<0.001) ){
    //   if(countnan>6){
    // 	checkmincurve=true;
    //   }
    //   countnan++;
    // }
     
    // if(RecoPar[3][size-1]<RecoPar[3][size-2] && RecoPar[3][size-2]<RecoPar[3][size-3] && checkmincurve==false){
    //   checkmincurve=false;
    // }
    iloop++;
  }
  
  int Nmin=0;
  int Nminsize=3;
  if(RecoPar[3].size()<3){
    Nminsize=RecoPar[3].size();
  }
  while(Nmin<Nminsize){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=RecoPar[2][FinalMinValueBin];
    GuessResultCor[Nmin][3]=RecoPar[3][FinalMinValueBin];
    GuessFinalMin[Nmin]=FinalMinValue;
    //cout<<"Guess Minimas 1st : Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }

  if(GuessResultCor[0][1]>180){
    GuessResultCor[0][1]=-360+GuessResultCor[0][1];
  }

  RecoPar[0].clear();
  RecoPar[1].clear();
  RecoPar[2].clear();
  RecoPar[3].clear();

  StartR=GuessResultCor[0][2]+250;
  StopR=GuessResultCor[0][2]-250;
  if(StartR<40){
    StartR=40;
  }
  if(StopR<40){
    StopR=40;
  }
  StepSizeR=50;
  checkmincurve=false;
  iloop=0;

  testTh=GuessResultCor[0][0];
  testPh=GuessResultCor[0][1];
  testR=StopR+1;    

  while(checkmincurve==false && testR>StopR){
   
    testTh=GuessResultCor[0][0];
    testPh=GuessResultCor[0][1];
    testR=StartR - iloop*StepSizeR ;
      
    Double_t ThPhR[3]={testTh*(Interferometer::pi/180),testPh*(Interferometer::pi/180),testR};
    Double_t XYZ[3]={0,0,0}; 
    Interferometer::ThPhRtoXYZ(ThPhR,XYZ);

    ParameterArray[iEnt]=XYZ[0];
    ParameterArray[iEnt+1]=XYZ[1];
    ParameterArray[iEnt+2]=XYZ[2];
    ParameterArray[iEnt+3]=testTh;
    ParameterArray[iEnt+4]=testPh;
    ParameterArray[iEnt+5]=testR;
    
    min=Interferometer::MinimizerThPh(testR, ParameterArray);  
      
    FinalTxCor[0]=ParameterArray[iEnt+6];
    FinalTxCor[1]=ParameterArray[iEnt+7];
    FinalTxCor[2]=testR;

    ThPhR[0]=FinalTxCor[0]*(Interferometer::pi/180.);
    ThPhR[1]=FinalTxCor[1]*(Interferometer::pi/180.);
    ThPhR[2]=FinalTxCor[2];

    if(std::isnan(min)==false && min!=1e9 && min<1e9 && FinalTxCor[0]>0.001 && min!=0){
      RecoPar[0].push_back(FinalTxCor[0]);
      RecoPar[1].push_back(FinalTxCor[1]);
      RecoPar[2].push_back(FinalTxCor[2]);
      RecoPar[3].push_back(min);      
    }

    int size=RecoPar[3].size();
    if(size>0 && iloop>=4){
      
      if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3] && iloop>=5) ){
    	//cout<<"found the minima!"<<endl;
    	checkmincurve=true;
      }

      if((std::isnan(min)==true || min==1e9 || min>1e9 || FinalTxCor[0]<0.001)){
	checkmincurve=true;
      }
    }
    
    // if(RecoPar[3][size-1]<RecoPar[3][size-2] && RecoPar[3][size-2]<RecoPar[3][size-3] && checkmincurve==false){
    // 	checkmincurve=false;
    // 	// if(iloop==7){
    // 	//   StepSizeR=StepSizeR*2;
    // 	// }	
    // }
    
    iloop++;
  }

  Nmin=0;
  Nminsize=3;
  if(RecoPar[3].size()<3){
    Nminsize=RecoPar[3].size();
  }
  while(Nmin<Nminsize){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=RecoPar[2][FinalMinValueBin];
    GuessResultCor[Nmin][3]=RecoPar[3][FinalMinValueBin];
    GuessFinalMin[Nmin]=FinalMinValue;
    //cout<<"Guess Minimas 2nd: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }

  if(GuessResultCor[0][1]>180){
    GuessResultCor[0][1]=-360+GuessResultCor[0][1];
  }

  RecoPar[0].clear();
  RecoPar[1].clear();
  RecoPar[2].clear();
  RecoPar[3].clear();

  min=0;
  
  StartR=4000;
  StopR=200;

  if(StopR<40){
    StopR=40;
  }
  
  StepSizeR=50;//(StopR-StartR)/20;
  checkmincurve=false;
  iloop=0;
  
  countnan=0;
  
  testTh=GuessResultCor[0][0];
  testPh=GuessResultCor[0][1];
  testR=StopR+1;
  
  while(checkmincurve==false && testR>StopR){
    
    testTh=GuessResultCor[0][0];
    testPh=GuessResultCor[0][1];
    testR=StartR - (iloop)*StepSizeR ;
      
    Double_t ThPhR[3]={testTh*(Interferometer::pi/180),testPh*(Interferometer::pi/180),testR};
    Double_t XYZ[3]={0,0,0}; 
    Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
    ParameterArray[iEnt]=XYZ[0];
    ParameterArray[iEnt+1]=XYZ[1];
    ParameterArray[iEnt+2]=XYZ[2];
    ParameterArray[iEnt+3]=testTh;
    ParameterArray[iEnt+4]=testPh;
    ParameterArray[iEnt+5]=testR;

    gsl_vector_set (ThPhRVec, 0, testTh);
    gsl_vector_set (ThPhRVec, 1, testPh);
    ParameterArray[iEnt+10]=testR;
      
    min=Interferometer::Minimizer_f(ThPhRVec, ParameterArray);

    FinalTxCor[0]=testTh;
    FinalTxCor[1]=testPh;
    FinalTxCor[2]=testR;
    
    if(std::isnan(min)==false && min!=1e9 && min<1e9 && FinalTxCor[0]>0.001 && min!=0){
      RecoPar[0].push_back(FinalTxCor[0]);
      RecoPar[1].push_back(FinalTxCor[1]);
      RecoPar[2].push_back(FinalTxCor[2]);
      RecoPar[3].push_back(min);
    }
    
    int size=RecoPar[3].size();
    if(size>0 && iloop>=5){
      //if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3]) ||  (fabs(RecoPar[3][size-1]-RecoPar[3][size-2])/RecoPar[3][size-1]<0.01 && fabs(RecoPar[3][size-3]-RecoPar[3][size-3])/RecoPar[3][size-2]<0.01) ){
      // if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3] && iloop>=5) ){
      //   cout<<"found the minima!"<<endl;
      //   checkmincurve=true;
      // }
    }
    // if((std::isnan(min)==true || min==1e9 || min>1e9 || FinalTxCor[0]<0.001) ){
    //   if(countnan>6){
    // 	checkmincurve=true;
    //   }
    //   countnan++;
    // }
     
    // if(RecoPar[3][size-1]<RecoPar[3][size-2] && RecoPar[3][size-2]<RecoPar[3][size-3] && checkmincurve==false){
    //   checkmincurve=false;
    // }
    iloop++;
  }
  
  Nmin=0;
  Nminsize=3;
  if(RecoPar[3].size()<3){
    Nminsize=RecoPar[3].size();
  }
  while(Nmin<Nminsize){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=RecoPar[2][FinalMinValueBin];
    GuessResultCor[Nmin][3]=RecoPar[3][FinalMinValueBin];
    GuessFinalMin[Nmin]=FinalMinValue;
    //cout<<"Guess Minimas 3rd : Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }

  if(GuessResultCor[0][1]>180){
    GuessResultCor[0][1]=-360+GuessResultCor[0][1];
  }

  RecoPar[0].clear();
  RecoPar[1].clear();
  RecoPar[2].clear();
  RecoPar[3].clear();

  StartR=GuessResultCor[0][2]+250;
  StopR=GuessResultCor[0][2]-250;
 
  if(StartR<40){
    StartR=40;
  }
  StepSizeR=100;
  checkmincurve=false;
  iloop=0;

  testTh=GuessResultCor[0][0];
  testPh=GuessResultCor[0][1];
  testR=StopR+1;    

  while(checkmincurve==false && testR>StopR){
   
    testTh=GuessResultCor[0][0];
    testPh=GuessResultCor[0][1];
    testR=StartR - iloop*StepSizeR ;
      
    Double_t ThPhR[3]={testTh*(Interferometer::pi/180),testPh*(Interferometer::pi/180),testR};
    Double_t XYZ[3]={0,0,0}; 
    Interferometer::ThPhRtoXYZ(ThPhR,XYZ);

    ParameterArray[iEnt]=XYZ[0];
    ParameterArray[iEnt+1]=XYZ[1];
    ParameterArray[iEnt+2]=XYZ[2];
    ParameterArray[iEnt+3]=testTh;
    ParameterArray[iEnt+4]=testPh;
    ParameterArray[iEnt+5]=testR;
    
    min=Interferometer::MinimizerThPh(testR, ParameterArray);  
      
    FinalTxCor[0]=ParameterArray[iEnt+6];
    FinalTxCor[1]=ParameterArray[iEnt+7];
    FinalTxCor[2]=testR;
  
    if(std::isnan(min)==false && min!=1e9 && min<1e9 && FinalTxCor[0]>0.001 &&  min!=0){
      RecoPar[0].push_back(FinalTxCor[0]);
      RecoPar[1].push_back(FinalTxCor[1]);
      RecoPar[2].push_back(FinalTxCor[2]);
      RecoPar[3].push_back(min);      
    }

    int size=RecoPar[3].size();
    if(size>0 && iloop>=4){
      //if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3]) ||  (fabs(RecoPar[3][size-1]-RecoPar[3][size-2])/RecoPar[3][size-1]<0.01 && fabs(RecoPar[3][size-3]-RecoPar[3][size-3])/RecoPar[3][size-2]<0.01) ){
      if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3] && iloop>=5) ){
    	//cout<<"found the minima!"<<endl;
    	checkmincurve=true;
      }

      if((std::isnan(min)==true || min==1e9 || min>1e9 || FinalTxCor[0]<0.001)){
  	checkmincurve=true;
      }
    }
    
    // if(RecoPar[3][size-1]<RecoPar[3][size-2] && RecoPar[3][size-2]<RecoPar[3][size-3] && checkmincurve==false){
    // 	checkmincurve=false;
    // 	// if(iloop==7){
    // 	//   StepSizeR=StepSizeR*2;
    // 	// }	
    // }
    
    iloop++;
  }

  Nmin=0;
  Nminsize=3;
  if(RecoPar[3].size()<3){
    Nminsize=RecoPar[3].size();
  }
  while(Nmin<Nminsize){
    int FinalMinValueBin=TMath::LocMin(RecoPar[3].size(),RecoPar[3].data());
    double FinalMinValue=RecoPar[3][FinalMinValueBin];	
    GuessResultCor[Nmin][0]=RecoPar[0][FinalMinValueBin];
    GuessResultCor[Nmin][1]=RecoPar[1][FinalMinValueBin];
    GuessResultCor[Nmin][2]=RecoPar[2][FinalMinValueBin];
    GuessResultCor[Nmin][3]=RecoPar[3][FinalMinValueBin];
    GuessFinalMin[Nmin]=FinalMinValue;
    //cout<<"Guess Minimas 4th: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }

  if(GuessResultCor[0][1]>180){
    GuessResultCor[0][1]=-360+GuessResultCor[0][1];
  }

  cout<<"Initial guess that will be passed to the minimizer: Theta "<<GuessResultCor[0][0]<<" ,Phi "<<GuessResultCor[0][1]<<" ,R "<<GuessResultCor[0][2]<<" , Fn value "<<GuessResultCor[0][3]<<endl;

  RecoPar[0].clear();
  RecoPar[1].clear();
  RecoPar[2].clear();
  RecoPar[3].clear();

}

void Interferometer::GetRecoFixedR(double InitialTxCor[3], double FinalTxCor[3], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2], vector <double> (&ChSNR)[2], double &FixedR, int IsItBelowStation, int max_iter){
  
  double ParameterArray[6*TotalAntennasRx+13];
  int iEnt=0;
  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChHitTime[iray][iRx];
      iEnt++;
    } 
  }  

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=IgnoreCh[iray][iRx];
      iEnt++;
    } 
  }

  for(int iray=0;iray<2;iray++){
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      ParameterArray[iEnt]=ChSNR[iray][iRx];
      iEnt++;
    }
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=0;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  ParameterArray[iEnt+13]=IsItBelowStation; 

  ParameterArray[iEnt+14]=max_iter;  

  double min=0;
  
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  double testTh=InitialTxCor[0];
  double testPh=InitialTxCor[1];
  double testR=FixedR;
      
  Double_t ThPhR[3]={testTh*(Interferometer::pi/180),testPh*(Interferometer::pi/180),testR};
  Double_t XYZ[3]={0,0,0}; 
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);

  if(testTh<45){
    ParameterArray[iEnt+12]=1;
  }
  
  ParameterArray[iEnt]=XYZ[0];
  ParameterArray[iEnt+1]=XYZ[1];
  ParameterArray[iEnt+2]=XYZ[2];
  ParameterArray[iEnt+3]=testTh;
  ParameterArray[iEnt+4]=testPh;
  ParameterArray[iEnt+5]=testR;
    
  min=Interferometer::MinimizerThPh(testR, ParameterArray);  
  
  FinalTxCor[0]=ParameterArray[iEnt+6];
  FinalTxCor[1]=ParameterArray[iEnt+7];
  FinalTxCor[2]=testR;
  
}

void Interferometer::DoInterferometery(double InitialTxCor[3], double FinalTxCor[3], vector<double> (&ChHitTime)[2], vector<int> (&IgnoreCh)[2],vector <double> (&ChSNR)[2], double &FinalMinValue, double &Duration,int &Iterations, double MinimizerRadialWidth,int IsItBelowStation, int max_iter){ 
  auto t1b = std::chrono::high_resolution_clock::now();
  
  Double_t ThPhR[3]={InitialTxCor[0],InitialTxCor[1],InitialTxCor[2]};
  Double_t XYZ[3]={0,0,0}; 
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
 
  ThPhR[0]=ThPhR[0]*(180./Interferometer::pi);
  ThPhR[1]=ThPhR[1]*(180./Interferometer::pi); 
  
  const int InitialValueNum=1;  
  vector <double> RecoXYZValues[3];
  vector <double> RecofMinValues;
  
  int FinalMinValueBin=0;
  double JitterNumber=5;  
    
  int iInVal=0;
  while(iInVal<InitialValueNum){
    double DummyRecoCor[3];
    double DummyMin=0;

    Interferometer::MinimizerThPhR(XYZ,ThPhR,DummyRecoCor, ChHitTime,IgnoreCh,ChSNR,DummyMin,Iterations,MinimizerRadialWidth,IsItBelowStation, max_iter);  
    
    RecoXYZValues[0].push_back(DummyRecoCor[0]);
    RecoXYZValues[1].push_back(DummyRecoCor[1]);
    RecoXYZValues[2].push_back(DummyRecoCor[2]);
    RecofMinValues.push_back(DummyMin);
    iInVal++;
  }
  
  if(RecofMinValues.size()!=0){
    FinalMinValueBin=TMath::LocMin(RecofMinValues.size(),RecofMinValues.data());
    FinalMinValue=RecofMinValues[FinalMinValueBin];
    FinalTxCor[0]=RecoXYZValues[0][FinalMinValueBin];
    FinalTxCor[1]=RecoXYZValues[1][FinalMinValueBin];
    FinalTxCor[2]=RecoXYZValues[2][FinalMinValueBin];
  }else{
    FinalTxCor[0]=0;
    FinalTxCor[1]=0;
    FinalTxCor[2]=0;
  }  

  if(FinalTxCor[1]>180){
    FinalTxCor[1]=360-FinalTxCor[1];
  }
  
  if(FinalTxCor[1]<-180){
    FinalTxCor[1]=360+FinalTxCor[1];
  }
  
  ThPhR[0]=FinalTxCor[0];
  ThPhR[1]=FinalTxCor[1];
  ThPhR[2]=FinalTxCor[2];
  
  auto t2b = std::chrono::high_resolution_clock::now();
  auto durationb = std::chrono::duration_cast<std::chrono::microseconds>( t2b - t1b ).count();
  Duration=durationb/1000;
}
