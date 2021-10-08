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

void Interferometer::GenerateChHitTimeAndCheckHits(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]){

  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]};
  
  int DHits=0,RHits=0;
  double timeD[TotalAntennasRx], timeR[TotalAntennasRx], timeRa[2][TotalAntennasRx];
  double RangD[TotalAntennasRx], RangR[TotalAntennasRx], RangRa[2][TotalAntennasRx];
  double RangRay[2][TotalAntennasRx];  
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing(0, TxDepth, Distance,RxDepth);
   
    // cout<<"*******For the Direct Ray********"<<endl;
    // cout<<"Launch Angle: "<<getresults[0]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[8]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[4]*pow(10,9)<<" ns"<<endl;
    // cout<<"*******For the Reflected Ray********"<<endl;
    // cout<<"Launch Angle: "<<getresults[1]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[9]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[5]*pow(10,9)<<" ns"<<endl;   
    // cout<<"Incident Angle in Ice on the Surface: "<<getresults[18]<<" deg"<<endl;
    // cout<<"*******For the Refracted[1] Ray********"<<endl;
    // cout<<"Launch Angle: "<<getresults[2]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[10]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[6]*pow(10,9)<<" ns"<<endl;
    // cout<<"*******For the Refracted[2] Ray********"<<endl;
    // cout<<"Launch Angle: "<<getresults[3]<<" deg"<<endl;
    // cout<<"Recieve Angle: "<<getresults[11]<<" deg"<<endl;
    // cout<<"Propogation Time: "<<getresults[7]*pow(10,9)<<" ns"<<endl;
    
    timeD[iRx]=RTresults[4]*pow(10,9);
    timeR[iRx]=RTresults[5]*pow(10,9);
    timeRa[0][iRx]=RTresults[6]*pow(10,9);
    timeRa[1][iRx]=RTresults[7]*pow(10,9);
    RangD[iRx]=RTresults[8];
    RangR[iRx]=RTresults[9];
    RangRa[0][iRx]=RTresults[10];
    RangRa[1][iRx]=RTresults[11];
    
    timeRay[0][iRx]=timeD[iRx];
    timeRay[1][iRx]=timeR[iRx];

    RangRay[0][iRx]=RangD[iRx];
    RangRay[1][iRx]=RangR[iRx];

    if(RangD[iRx]!=-1000){
      timeRay[0][iRx]=timeD[iRx];
      RangRay[0][iRx]=RangD[iRx];
    }

    if(RangR[iRx]!=-1000){
      timeRay[1][iRx]=timeR[iRx];
      RangRay[1][iRx]=RangR[iRx];
    }
    
    if(RangRa[0][iRx]!=-1000 && RangD[iRx]!=-1000){
      timeRay[1][iRx]=timeRa[0][iRx];
      RangRay[1][iRx]=RangRa[0][iRx];
    }

    if(RangRa[0][iRx]!=-1000 && RangR[iRx]!=-1000){
      timeRay[0][iRx]=timeRa[0][iRx];
      RangRay[0][iRx]=RangRa[0][iRx];
    }
    
    if(RangRa[0][iRx]!=-1000 && RangD[iRx]==-1000 && RangR[iRx]==-1000){
      timeRay[0][iRx]=timeRa[0][iRx];
      RangRay[0][iRx]=RangRa[0][iRx];
    }

    if(RangRa[1][iRx]!=-1000 && RangD[iRx]==-1000 && RangR[iRx]==-1000){
      timeRay[1][iRx]=timeRa[1][iRx];
      RangRay[1][iRx]=RangRa[1][iRx];
    }

    if(timeRay[0][iRx]>timeRay[1][iRx]){
      swap(timeRay[0][iRx],timeRay[1][iRx]);
      swap(RangRay[0][iRx],RangRay[1][iRx]);
    }
    
    if(RangRay[0][iRx]==-1000){
      IgnoreCh[0][iRx]=0;
    }

    if(RangRay[1][iRx]==-1000){
      IgnoreCh[1][iRx]=0;
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

void Interferometer::GenerateChHitTimeAndCheckHits_Cnz(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]){

  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]};
  
  int DHits=0,RHits=0;
  // double LangD=0;double LangR=0;double LangRa=0;
  double timeD[TotalAntennasRx], timeR[TotalAntennasRx], timeRa[TotalAntennasRx];
  double RangD[TotalAntennasRx], RangR[TotalAntennasRx], RangRa[TotalAntennasRx];
  double RangRay[2][TotalAntennasRx];  
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxDepth=AntennaCoordTx[2];
    double *RTresults=IceRayTracing::IceRayTracing_Cnz(0, TxDepth, Distance,RxDepth,(IceRayTracing::Getnz(TxDepth)+IceRayTracing::Getnz(RxDepth))/2);

    timeD[iRx]=RTresults[2]*pow(10,9);
    timeR[iRx]=RTresults[3]*pow(10,9);
   
    RangD[iRx]=RTresults[4];
    RangR[iRx]=RTresults[5];
    
    timeRay[0][iRx]=timeD[iRx];
    timeRay[1][iRx]=timeR[iRx];

    RangRay[0][iRx]=RangD[iRx];
    RangRay[1][iRx]=RangR[iRx];
    
    DHits++;
    RHits++;
    
    delete [] RTresults;  
  }  
  //cout<<"hits are "<<DHits<<" "<<RHits<<endl; 
}


void Interferometer::GenerateChHitTimeAndCheckHits_Air(double TxCor[3],double timeRay[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx]){

  double AntennaCoordTx[3]={TxCor[0],TxCor[1],TxCor[2]};
  
  int DHits=0;
  // double LangD=0;;
  double timeD[TotalAntennasRx];
  double RangD[TotalAntennasRx];
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    double Distance=sqrt(pow(AntennaCoordRx[iRx][0]-AntennaCoordTx[0],2) + pow(AntennaCoordRx[iRx][1]-AntennaCoordTx[1],2));
    double RxDepth=AntennaCoordRx[iRx][2];
    double TxHeight=AntennaCoordTx[2];   
    //cout<<"RT par are "<<RxDepth<<" "<<Distance<<" "<<TxHeight<<endl;
    double *RTresults=IceRayTracing::GetDirectRayPar_Air(RxDepth, Distance, TxHeight);
    timeD[iRx]=RTresults[2]*pow(10,9);
    RangD[iRx]=RTresults[0];
    timeRay[0][iRx]=timeD[iRx];

    if(RangD[iRx]==-1000){
      IgnoreCh[0][iRx]=0;
    }

    IgnoreCh[1][iRx]=0;
    
    DHits++;   
    
    delete [] RTresults;  
  }
}

bool Interferometer::CheckTrigger(int IgnoreCh[2][TotalAntennasRx]){

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
    if(count[iray][0]+count[iray][1]>=3){
      if(count[iray][0]>=4){
	RayTrigger[iray]=true;
      } 
      if(count[iray][1]>=4){
	RayTrigger[iray]=true;
      }
    }
  }

  bool DidStationTrigger=false;

  if(RayTrigger[0]==true || RayTrigger[1]==true){
    DidStationTrigger=true;
  }

  if(RayTrigger[0]==false && RayTrigger[1]==false){
    DidStationTrigger=false;
  }

  return DidStationTrigger;
}

void Interferometer::ReadChHitTimeFromData(const char * filename,double ChHitTime[2][TotalAntennasRx]){

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

void Interferometer::FindFirstHitAndNormalizeHitTime(double ChHitTime[2][TotalAntennasRx],int IgnoreCh[2][TotalAntennasRx],double ChDRTime[TotalAntennasRx], int ChHitOrder[TotalAntennasRx]){

  double FirstHitTime[2]={0,0};
  int FirstHitCh[2]={0,0};

  double LastHitTime[2]={0,0};
  int LastHitCh[2]={0,0};
  vector <double> ChHitTimeB[2];
  vector <double> ChHitTimeC[2];
  vector <int> ChHitNum[2];

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ChHitOrder[iRx]=-1;
  }
  
  for(int iray=0;iray<2;iray++){ 
    for(int iRx=0;iRx<TotalAntennasRx;iRx++){
      if(IgnoreCh[iray][iRx]==1){
	ChHitTimeB[iray].push_back(ChHitTime[iray][iRx]);
	ChHitNum[iray].push_back(iRx);
      }
    }
  }
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(IgnoreCh[0][iRx]==1 && IgnoreCh[1][iRx]==1){
      ChDRTime[iRx]=ChHitTime[1][iRx]-ChHitTime[0][iRx];
    }else{
      ChDRTime[iRx]=0;
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
	  ChHitTime[iray][iRx]=ChHitTime[iray][iRx]-FirstHitTime[0];
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

  if(ChHitTimeB[0].size()!=0){
    for(int iRx=0;iRx<ChHitTimeB[0].size();iRx++){
      int HitBin=TMath::LocMin(ChHitTimeB[0].size(),ChHitTimeB[0].data());
      int HitCh=ChHitNum[0][HitBin];
      ChHitOrder[iRx]=HitCh;
      ChHitTimeB[0][HitBin]=1e9;
    }
  }  
  
  // for(int iRx=0;iRx<TotalAntennasRx;iRx++){
  //   //ChHitTime[iray][iRx]=0;
  //   cout<<iRx<<" "<<ChHitTime[0][iRx]<<" "<<ChHitTime[0][iRx]<<endl;
  // }
}

void Interferometer::AddGaussianJitterToHitTimes(double JitterNumber,double ChHitTime[2][TotalAntennasRx]){

  ////Introduce jitter +-JitterNumber
  TRandom3 *RandNumIni = new TRandom3(0);
  for (int iRx=0; iRx<TotalAntennasRx; iRx++){
    for(int iray=0;iray<2;iray++){ 
      //double RandNum = (RandNumIni->Rndm(iRx)*2*JitterNumber)-JitterNumber;
      double RandNum = (RandNumIni->Gaus(0,0.4));
      if(RandNum>=JitterNumber){
      	RandNum=JitterNumber;
      }
      if(RandNum<=-JitterNumber){
      	RandNum=-JitterNumber;
      }
      ChHitTime[iray][iRx]=ChHitTime[iray][iRx]+RandNum;
    }
  }
  delete RandNumIni;
}

double Interferometer::Minimizer_f(const gsl_vector *v, void *params){

  double *p = (double *)params;
  
  double ExpectedX=p[7*TotalAntennasRx];
  double ExpectedY=p[7*TotalAntennasRx+1];
  double ExpectedZ=p[7*TotalAntennasRx+2];

  double ExpectedTh=p[7*TotalAntennasRx+3];
  double ExpectedPh=p[7*TotalAntennasRx+4];
  double ExpectedR=p[7*TotalAntennasRx+5];
  
  double theta, phi, r;
  theta = gsl_vector_get(v, 0)*(Interferometer::pi/180.0);
  phi = gsl_vector_get(v, 1)*(Interferometer::pi/180.0);
  r = p[7*TotalAntennasRx+10];
  
  Double_t ThPhR[3]={theta,phi,r};
  Double_t XYZ[3]={0,0,0};
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
  
  double output=0;
  double chi2=0,chi2d=0;

  double AntennaCoordTx[3]={0,0,0};
  for(int ixyz=0;ixyz<3;ixyz++){
    AntennaCoordTx[ixyz]=XYZ[ixyz];
  }
    
  double timeRay[2][TotalAntennasRx];
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
    
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
    if(p[2*TotalAntennasRx +iRx]==0){
      IgnoreCh[0][iRx]=0;
    }
    if(p[3*TotalAntennasRx +iRx]==0){
      IgnoreCh[1][iRx]=0;
    }
    ChDRTime[iRx]=0;
  }    

  if(XYZ[2]<0 || p[7*TotalAntennasRx+12]==1){
    // cout<<"working for ice "<<endl;
    Interferometer::GenerateChHitTimeAndCheckHits(AntennaCoordTx,timeRay,IgnoreCh);  
  }else{
    //cout<<"working for air "<<endl;
    Interferometer::GenerateChHitTimeAndCheckHits_Air(AntennaCoordTx,timeRay,IgnoreCh);  
  }

  Interferometer::FindFirstHitAndNormalizeHitTime(timeRay,IgnoreCh,ChDRTime,ChHitOrder);
  
  Double_t SumSNRD=0, SumSNRR=0,SumSNR=0;
  int chanD[2]={0,0};
  int chanR[2]={0,0};
  int chanDsame=0,chanRsame=0;
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
  }
  SumSNR=SumSNRD+SumSNRR;
 
  // double ChHitOrderDiff=0;
  // double ChHitOrderDiffTot=0;
  // for(int iRx=0;iRx<TotalAntennasRx;iRx++){
  //   if(p[6*TotalAntennasRx +iRx]!=-1  && ChHitOrder[iRx]!=-1){ 
  // 	if(p[6*TotalAntennasRx +iRx]!=ChHitOrder[iRx]){
  // 	  ChHitOrderDiff+=100;
  // 	}
  // 	ChHitOrderDiffTot+=1000;
  //   }
  // }
  //ChHitOrderDiff=ChHitOrderDiff/ChHitOrderDiffTot;

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){ 
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
      chi2+=pow((p[4*TotalAntennasRx +iRx]*(timeRay[0][iRx] - p[0+iRx]))/SumSNR,2);
      chi2d+=pow(p[4*TotalAntennasRx +iRx]/SumSNR,2);
    }
    if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      chi2+=pow((p[5*TotalAntennasRx +iRx]*(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx]))/SumSNR,2);
      chi2d+=pow(p[5*TotalAntennasRx +iRx]/SumSNR,2);
    }
  }
   
  
  if(chanD[1]>=chanDsame && chanR[1]>=chanRsame && chanDsame>=chanD[0] && chanRsame>=chanR[0]){   
    output=(chi2/chi2d);
    if(std::isnan(output)){
      output=1e9;
    }
  }else{
    //cout<<chanD[1]<<" "<<chanDsame<<" "<<chanR[1]<<" "<<chanRsame<<" "<<chanDsame<<" "<<chanD[0]<<" "<<chanRsame<<" "<<chanR[0]<<endl;
    output=1e9+pow(chi2/chi2d,2); 
  }
  
  if((gsl_vector_get(v, 0)>179.9 || gsl_vector_get(v, 0)<0.1)  && p[7*TotalAntennasRx+12]==0){
    output=1e9+pow(chi2/chi2d,2);
  }

  if((gsl_vector_get(v, 0)<90.2 || gsl_vector_get(v, 0)>179.9 || gsl_vector_get(v, 0)<0.1) && p[7*TotalAntennasRx+12]==1){
    output=1e9+pow(chi2/chi2d,2);
  }

  if(std::isnan(output)){
    output=1e9;
  }
  
  return output;
}

double Interferometer::Minimizer_fCnz(const gsl_vector *v, void *params){

  double *p = (double *)params; 

  double ExpectedX=p[7*TotalAntennasRx];
  double ExpectedY=p[7*TotalAntennasRx+1];
  double ExpectedZ=p[7*TotalAntennasRx+2];
  double distance=p[7*TotalAntennasRx+5];

  double theta, phi, r;
  theta = gsl_vector_get(v, 0)*(Interferometer::pi/180.0);
  phi = gsl_vector_get(v, 1)*(Interferometer::pi/180.0);
  r = p[7*TotalAntennasRx+10];

  Double_t ThPhR[3]={theta,phi,r};
  Double_t XYZ[3]={0,0,0};
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
  
  double output=0;
  double chi2=0,chi2d=0;
  
  double AntennaCoordTx[3]={0,0,0};
  for(int ixyz=0;ixyz<3;ixyz++){
    AntennaCoordTx[ixyz]=XYZ[ixyz];
  }
    
  double timeRay[2][TotalAntennasRx];
  int IgnoreCh[2][TotalAntennasRx];
  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    for(int iray=0;iray<2;iray++){
      IgnoreCh[iray][iRx]=1;
    }
    if(p[2*TotalAntennasRx +iRx]==0){
      IgnoreCh[0][iRx]=0;
    }
    if(p[3*TotalAntennasRx +iRx]==0){
      IgnoreCh[1][iRx]=0;
    }

    ChDRTime[iRx]=0;
  }

  if(XYZ[2]<0 || p[7*TotalAntennasRx+12]==1){
    Interferometer::GenerateChHitTimeAndCheckHits_Cnz(AntennaCoordTx,timeRay,IgnoreCh);  
  }else{
    Interferometer::GenerateChHitTimeAndCheckHits_Air(AntennaCoordTx,timeRay,IgnoreCh);  
  }
    
  Interferometer::FindFirstHitAndNormalizeHitTime(timeRay,IgnoreCh,ChDRTime,ChHitOrder);

  Double_t SumSNRD=0, SumSNRR=0,SumSNR=0;
  int chanD[2]={0,0};
  int chanR[2]={0,0};
  int chanDsame=0,chanRsame=0;
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
  }
  SumSNR=SumSNRD+SumSNRR;

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(p[2*TotalAntennasRx +iRx]!=0 && IgnoreCh[0][iRx]!=0){
      chi2+=pow(((p[4*TotalAntennasRx +iRx]*(timeRay[0][iRx] - p[0+iRx]))/SumSNR),2);
      chi2d+=pow(p[4*TotalAntennasRx +iRx]/SumSNR,2);
	
    }
    if(p[3*TotalAntennasRx +iRx]!=0 && IgnoreCh[1][iRx]!=0){
      chi2+=pow(((p[5*TotalAntennasRx +iRx]*(timeRay[1][iRx] - p[1*TotalAntennasRx+iRx]))/SumSNR),2);
      chi2d+=pow(p[5*TotalAntennasRx +iRx]/SumSNR,2);
    }
  }
  
  if(chanD[1]>=chanDsame && chanR[1]>=chanRsame && chanDsame>=chanD[0] && chanRsame>=chanR[0]){
    output=(chi2/chi2d);
    if(std::isnan(output)){
      output=1e9;
    }   
  }else{
    output=1e9+pow(chi2/chi2d,2); 
  }
  
  if((gsl_vector_get(v, 0)>179.9 || gsl_vector_get(v, 0)<0.1)  && p[7*TotalAntennasRx+12]==0){
    output=1e9+pow(chi2/chi2d,2);
  }

  if((gsl_vector_get(v, 0)<90.2 || gsl_vector_get(v, 0)>179.9 || gsl_vector_get(v, 0)<0.1) && p[7*TotalAntennasRx+12]==1){
    output=1e9+pow(chi2/chi2d,2);
  }
  
  if(std::isnan(output)){
    output=1e9;
  }  

  return output;
}

double Interferometer::MinimizerThPh(double x, void * params)
{
 
  double *p = (double *)params;
  
  p[7*TotalAntennasRx+10]=x;
  
  const gsl_multimin_fminimizer_type *MinimisationType = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *MinimizerWorkSpace = NULL;
  gsl_vector *XYZStepSizeVector, *XYZVec;
  gsl_multimin_function MultiDimMinFunc;

  size_t Iter = 0;
  int Status;
  double Size;
    
  /* Starting point */
  XYZVec = gsl_vector_alloc (2);
  gsl_vector_set (XYZVec, 0, p[7*TotalAntennasRx+3]);
  gsl_vector_set (XYZVec, 1, p[7*TotalAntennasRx+4]);
  
  double ExpectedDisplacement=p[7*TotalAntennasRx+5];
  double ExpectedUncertainty=p[7*TotalAntennasRx+9];
  double ThetaPhiStepSize=atan(ExpectedUncertainty/ExpectedDisplacement)*(180./Interferometer::pi);
  
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

  double FinalMinValue=0;
  double FinalTxCor[2]; 
  int Iterations=0;
  int Max_Iter=50;
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
      // printf ("%5d %7.3f %7.3f f() = %7.3f size = %.3f\n",
      //         Iter,
      //         gsl_vector_get (MinimizerWorkSpace->x, 0),
      //         gsl_vector_get (MinimizerWorkSpace->x, 1),
      //         MinimizerWorkSpace->fval,Size);
     
      //cout<<Iter<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 0)<<" "<<gsl_vector_get (MinimizerWorkSpace->x, 1) <<" "<<MinimizerWorkSpace->fval<<" "<<Size<<" "<<(*((MultiDimMinFunc).f))(XYZVec,(MultiDimMinFunc).params)<<endl;
      
      FinalMinValue=MinimizerWorkSpace->fval;
      FinalTxCor[0]=gsl_vector_get (MinimizerWorkSpace->x, 0);
      FinalTxCor[1]=gsl_vector_get (MinimizerWorkSpace->x, 1);
    }
  while (Status == GSL_CONTINUE && Iter < Max_Iter);

  //printf ("error: %s\n", gsl_strerror (Status));
  p[7*TotalAntennasRx+11]=FinalMinValue;
  
  p[7*TotalAntennasRx+6]=FinalTxCor[0];
  p[7*TotalAntennasRx+7]=FinalTxCor[1];
  p[7*TotalAntennasRx+8]=x;
  
  if(FinalMinValue==0 ){
    p[7*TotalAntennasRx+11]=1e9;
    FinalMinValue=1e9;
  }
  
  gsl_vector_free(XYZVec);
  gsl_vector_free(XYZStepSizeVector);
  gsl_multimin_fminimizer_free (MinimizerWorkSpace);
  
  return FinalMinValue;
}

void Interferometer::MinimizerThPhR(double InitialTxCor_XYZ[3], double InitialTxCor_ThPhR[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int ChHitOrder[TotalAntennasRx], double &FinalMinValue, int &Iterations, double MinimizerRadialWidth){
 
  double ParameterArray[7*TotalAntennasRx+13];
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

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
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

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;

  if(InitialTxCor_ThPhR[0]>90){
    ParameterArray[iEnt+12]=1;
  }
  
  int status;
  int iter = 0, max_iter = 50;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  double m = InitialTxCor_ThPhR[2];
  double a = InitialTxCor_ThPhR[2]-MinimizerRadialWidth, b = InitialTxCor_ThPhR[2]+MinimizerRadialWidth;
  if(a<0){
    a=10;
  }
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
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_min_fminimizer_free (s);
  FinalTxCor[0]=ParameterArray[iEnt+6];
  FinalTxCor[1]=ParameterArray[iEnt+7];
  FinalTxCor[2]=m;
  FinalMinValue=ParameterArray[7*TotalAntennasRx+11];
  
}

void Interferometer::RootThPhR(double InitialTxCor_XYZ[3], double InitialTxCor_ThPhR[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], int ChHitOrder[TotalAntennasRx], double &FinalMinValue, int &Iterations){
  
  double ParameterArray[7*TotalAntennasRx+13];
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

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
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

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  
  int status1,status2;
  int iter = 0, max_iter = 10;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = InitialTxCor_ThPhR[2];
  double x_lo = InitialTxCor_ThPhR[2]-100, x_hi = InitialTxCor_ThPhR[2]+100;
  if(x_lo<0){
    x_lo=10;
  }
 
  gsl_function F;
  
  F.function = &MinimizerThPh;
  F.params = ParameterArray;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  // printf ("using %s method\n",
  //         gsl_root_fsolver_name (s));

  // printf ("%5s [%9s, %9s] %9s %10s %9s\n",
  //         "iter", "lower", "upper", "root",
  //         "err", "err(est)");

  do
    {
      iter++;
      status1 = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      
      double checkval=(*((F).function))(r,(F).params);  
      status1=gsl_root_test_residual(checkval, 2);
      status2= gsl_root_test_interval (x_lo, x_hi,0.01, 0);
      
      // if (status1 == GSL_SUCCESS && status2 == GSL_SUCCESS)
      //   printf ("Converged:\n");

      // printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
      //         iter, x_lo, x_hi,
      //         r, r,
      //         x_hi - x_lo);
    }
  while ((status1 == GSL_CONTINUE || status2 == GSL_CONTINUE) && iter < max_iter);

  gsl_root_fsolver_free (s);

  FinalTxCor[0]=ParameterArray[iEnt+6];
  FinalTxCor[1]=ParameterArray[iEnt+7];
  FinalTxCor[2]=r;
  FinalMinValue=ParameterArray[7*TotalAntennasRx+11];
 
}

double Interferometer::GetChiSquaredThPhR(double UserCor[3] , double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx]){

  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
  
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime,ChHitOrder);

  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  double ParameterArray[7*TotalAntennasRx+13];
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
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;

  double Tht=UserCor[0];
  double Pht=UserCor[1];
  double Rt=UserCor[2];
    
  if(Tht>90){
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
  min=Interferometer::Minimizer_f(ThPhRVec, ParameterArray);
    
  gsl_vector_free(ThPhRVec);

  return min;
  
}


void Interferometer::SearchApproxiMin(int C_nz, double StartCor[3],double GuessResultCor[3][3],double ParameterArray[7*TotalAntennasRx+13],int &iEnt, double StartDistance, bool CheckBelowSurface){

  double Number=0;
  if(CheckBelowSurface==true){
    Number=90;
  }
  Double_t NumBinsTh=20,NumBinsPh=20,NumBinsR=4;
  Double_t StartTh=2+Number,StartPh=-180,StartR=StartDistance+20;
  Double_t StopTh=178,StopPh=179,StopR=StartDistance+100;  
  if(C_nz==0){
    NumBinsTh=5,NumBinsPh=5,NumBinsR=4;
    StartTh=StartCor[0]-20,StartPh=StartCor[1]-20,StartR=StartDistance+20;
    StopTh=StartCor[0]+20,StopPh=StartCor[1]+20,StopR=StartDistance+100;
  }
  if(StartTh<=0){
    StartTh=2;
  } 
  if(StopTh>=180){
    StopTh=178;
  }
  if(StopTh>178 && CheckBelowSurface==true){
    StopTh=178;
  }
  
  if(StartR<0){
    StartR=10;
  }
 
  Double_t StepSizeTh=(StopTh-StartTh)/NumBinsTh,StepSizePh=(StopPh-StartPh)/NumBinsPh,StepSizeR=(StopR-StartR)/NumBinsR;
  StepSizeR=20;
  NumBinsR=((StopR-StartR)/StepSizeR)+1;
  vector <double> RecoPar[4];
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  for(double i=0; i<=NumBinsTh;i++){
    double Tht=i*StepSizeTh + StartTh;
    for(double j=0; j<=NumBinsPh;j++){	
      double Pht=j*StepSizePh + StartPh;
      for(double k=0; k<=NumBinsR;k++){
  	double Rt=k*StepSizeR + StartR;

	if(Tht==90 && CheckBelowSurface==true){
	  Tht=90.2;
	}

	if(Tht==90 && CheckBelowSurface==false){
	  Tht=89.8;
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
  	if(C_nz==0){
  	  min=Interferometer::Minimizer_f(ThPhRVec, ParameterArray);
  	}else{
  	  min=Interferometer::Minimizer_fCnz(ThPhRVec, ParameterArray);
  	}
	
  	if(std::isnan(min)==false && min!=1e9 && min<1e9 && Tht!=90){
  	  //cout<<Tht<<" "<<Pht<<" "<<Rt<<" "<<min<<endl;
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
    cout<<"Guess Minimas: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }

  if(C_nz==1){
    StartCor[0]=GuessResultCor[0][0];
    StartCor[1]=GuessResultCor[0][1];
    StartCor[2]=GuessResultCor[0][2];
  }else{
    StartCor[0]=0;
    StartCor[1]=0;
    StartCor[2]=0;
  }
  
}

void Interferometer::GetApproximateMinUserCor(vector <double> UserCor[3] ,double GuessResultCor[3][3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx]){

  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
  
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime,ChHitOrder);

  vector <double> RecoPar[4];
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  double ParameterArray[7*TotalAntennasRx+13];
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
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  
  for(int i=0;i<UserCor[0].size();i++){
    double Tht=UserCor[0][i];
    double Pht=UserCor[1][i];
    double Rt=UserCor[2][i];
    
    if(Tht>90){
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
    cout<<"Guess Minimas: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }
}

void Interferometer::GetApproximateMinThPhR(double GuessResultCor[3][3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], double StartDistance){
    
  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
  
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime,ChHitOrder);
  
  double ParameterArray[7*TotalAntennasRx+13];
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
  
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  
  bool CheckBelowSurface=false;
  bool IsItBelowStation=true;
  
  double StartCor[3];
  Interferometer::SearchApproxiMin(1,StartCor,GuessResultCor,ParameterArray,iEnt,StartDistance,CheckBelowSurface);

  double XYZ[3]={0,0,0};
  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    if(ChHitTime[0][iRx]==0 && IgnoreCh[0][iRx]==1){
      XYZ[0]=AntennaCoordRx[iRx][0]-AvgAntennaCoordRx[0];
      XYZ[1]=AntennaCoordRx[iRx][1]-AvgAntennaCoordRx[1];
      XYZ[2]=AntennaCoordRx[iRx][2]-AvgAntennaCoordRx[2];

      if(XYZ[2]>0){
	cout<<"Ch."<<iRx<<" was hit first and it is definitely above the station"<<endl;
	IsItBelowStation=false;
      }
      if(XYZ[2]<0){
	cout<<"Ch."<<iRx<<" was hit first and it is definitely below the station in ice"<<endl;
	IsItBelowStation=true;
      }
    }
  }
  
  if(IsItBelowStation==false){

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
      cout<<"It could be in ice or in air"<<endl;
      ParameterArray[iEnt+12]=0;
    }

    if(AreAllRIgnored==false && IsdtDRtimeInWindow==true){
      CheckBelowSurface=true;
      cout<<"It is definitely in ice"<<endl;
      ParameterArray[iEnt+12]=1;
      Interferometer::SearchApproxiMin(1,StartCor,GuessResultCor,ParameterArray,iEnt,StartDistance,CheckBelowSurface);
    }
    
  }

  if(IsItBelowStation==true){
    CheckBelowSurface=true;
    ParameterArray[iEnt+12]=1;
    Interferometer::SearchApproxiMin(1,StartCor,GuessResultCor,ParameterArray,iEnt,StartDistance,CheckBelowSurface);
  }
  
  Interferometer::SearchApproxiMin(0,StartCor,GuessResultCor,ParameterArray,iEnt,StartDistance,CheckBelowSurface);

}

void Interferometer::GetApproximateDistance(double GuessResultCor[3][3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx]){
    
  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
  vector <double> RecoPar[4];
  
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime,ChHitOrder);
  
  double ParameterArray[7*TotalAntennasRx+13];
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

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;

  if(GuessResultCor[0][0]>90){
    ParameterArray[iEnt+12]=1;
  }
  double StartR=GuessResultCor[0][2];
  double StopR=500+GuessResultCor[0][2];
  double StepSizeR=100;//(StopR-StartR)/20;
  double FinalTxCor[3];
  bool checkmincurve=false;
  int iloop=0;
  double min=0;
  
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  while(checkmincurve==false){
    
    if(RecoPar[3].size()==0){      
      double testTh=GuessResultCor[0][0];
      double testPh=GuessResultCor[0][1];
      double testR=iloop*StepSizeR + StartR;
      
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
    }
    
    if(iloop>0 && RecoPar[3].size()>0){
      double testTh=RecoPar[0][0];
      double testPh=RecoPar[1][0];
      double testR=(iloop)*StepSizeR + StartR;
      
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
    }

    if(std::isnan(min)==false && min!=1e9 && min<1e9 && FinalTxCor[0]>0.001){
      RecoPar[0].push_back(FinalTxCor[0]);
      RecoPar[1].push_back(FinalTxCor[1]);
      RecoPar[2].push_back(FinalTxCor[2]);
      RecoPar[3].push_back(min);
    }
    
    int size=RecoPar[3].size();
    if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3] && iloop>=5) ||  (fabs(RecoPar[3][size-1]-RecoPar[3][size-2])/RecoPar[3][size-1]<0.01 && iloop>=5) ){
      //cout<<"found the minima!"<<endl;
      checkmincurve=true;
    }
     
    if((std::isnan(min)==true || min==1e9 || min>1e9 || FinalTxCor[0]<0.001) && iloop>=5 ){
      checkmincurve=true;
    }
     
    if(RecoPar[3][size-1]<RecoPar[3][size-2] && RecoPar[3][size-2]<RecoPar[3][size-3] && checkmincurve==false){
      checkmincurve=false;
    }
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
    cout<<"Guess Minimas 1st : Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
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
  
  StartR=10;
  StopR=GuessResultCor[0][2];
  StepSizeR=(StopR-StartR)/5;
  checkmincurve=false;
  iloop=0;
  min=0;

  while(checkmincurve==false){
      
    double testTh=GuessResultCor[0][0];
    double testPh=GuessResultCor[0][1];
    double testR=iloop*StepSizeR + StartR;
      
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
    
    if(std::isnan(min)==false && min!=1e9 && min<1e9 && FinalTxCor[0]>0.001){
      RecoPar[0].push_back(FinalTxCor[0]);
      RecoPar[1].push_back(FinalTxCor[1]);
      RecoPar[2].push_back(FinalTxCor[2]);
      RecoPar[3].push_back(min);      
    }

    int size=RecoPar[3].size();
    if((RecoPar[3][size-1]>RecoPar[3][size-2] && RecoPar[3][size-2]>RecoPar[3][size-3] && iloop>=5) ||  (fabs(RecoPar[3][size-1]-RecoPar[3][size-2])/RecoPar[3][size-1]<0.01 && iloop>=5) ){
      //cout<<"found the minima!"<<endl;
      checkmincurve=true;
    }

    if((std::isnan(min)==true || min==1e9 || min>1e9 || FinalTxCor[0]<0.001) && iloop>=5 ){
      checkmincurve=true;
    }
    
    if(RecoPar[3][size-1]<RecoPar[3][size-2] && RecoPar[3][size-2]<RecoPar[3][size-3] && checkmincurve==false){
      checkmincurve=false;
      if(iloop==7){
	StepSizeR=StepSizeR*2;
      }	
    }
    
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
    cout<<"Guess Minimas 2nd: Tht "<<GuessResultCor[Nmin][0]<<" Pht "<<GuessResultCor[Nmin][1]<<" Rt "<<GuessResultCor[Nmin][2]<<" min "<<FinalMinValue<<endl;
    RecoPar[3][FinalMinValueBin]=10000000000;    
    Nmin++;
  }
 
  if(GuessResultCor[0][1]>180){
    GuessResultCor[0][1]=-360+GuessResultCor[0][1];
  }

  if(fabs(GuessResultCor[0][0]-90)<0.001 && ParameterArray[iEnt+12]==1){
    GuessResultCor[0][0]=90.5;
  }
  
 
}

void Interferometer::GetRecoFixedR(double InitialTxCor[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx], double ChSNR[2][TotalAntennasRx], double &FixedR){
    
  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
  
  //Interferometer::AddGaussianJitterToHitTimes(10,ChHitTime);
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime,ChHitOrder);
  
  double ParameterArray[7*TotalAntennasRx+13];
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

  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
    ParameterArray[iEnt]=ChHitOrder[iRx];
    iEnt++;
  }
  
  ParameterArray[iEnt+6]=0;
  ParameterArray[iEnt+7]=0;
  ParameterArray[iEnt+8]=0;

  ParameterArray[iEnt+9]=ExpectedUncertainty;
  ParameterArray[iEnt+10]=0;
  ParameterArray[iEnt+11]=0;

  ParameterArray[iEnt+12]=0;
  
  double min=0;
  
  gsl_vector *ThPhRVec;
  ThPhRVec = gsl_vector_alloc (2);
  
  double testTh=InitialTxCor[0];
  double testPh=InitialTxCor[1];
  double testR=FixedR;
      
  Double_t ThPhR[3]={testTh*(Interferometer::pi/180),testPh*(Interferometer::pi/180),testR};
  Double_t XYZ[3]={0,0,0}; 
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);

  if(testTh>90){
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


void Interferometer::DoInterferometery(double InitialTxCor[3], double FinalTxCor[3], double ExpectedUncertainty, double ChHitTime[2][TotalAntennasRx], int IgnoreCh[2][TotalAntennasRx],double ChSNR[2][TotalAntennasRx], double &FinalMinValue, double &Duration,int &Iterations, double MinimizerRadialWidth){ 
  auto t1b = std::chrono::high_resolution_clock::now();
  
  Double_t ThPhR[3]={InitialTxCor[0],InitialTxCor[1],InitialTxCor[2]};
  Double_t XYZ[3]={0,0,0}; 
  Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
  ThPhR[0]=(ThPhR[0])*(180./Interferometer::pi);
  ThPhR[1]=ThPhR[1]*(180./Interferometer::pi); 
  
  const int InitialValueNum=1;  
  vector <double> RecoXYZValues[3];
  vector <double> RecofMinValues;
  
  int FinalMinValueBin=0;
  double JitterNumber=5; 
  double ChDRTime[TotalAntennasRx];
  int ChHitOrder[TotalAntennasRx];
  //Interferometer::AddGaussianJitterToHitTimes(JitterNumber,ChHitTime);
  Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime,ChHitOrder);
    
  int iInVal=0;
  while(iInVal<InitialValueNum){
    double DummyRecoCor[3];
    double DummyMin=0;

    //Interferometer::RootThPhR(XYZ,ThPhR,DummyRecoCor,ExpectedUncertainty,ChHitTime,IgnoreCh,ChSNR,ChHitOrder,DummyMin,Iterations);
    Interferometer::MinimizerThPhR(XYZ,ThPhR,DummyRecoCor,ExpectedUncertainty,ChHitTime,IgnoreCh,ChSNR,ChHitOrder,DummyMin,Iterations,MinimizerRadialWidth);  

    // //if(DummyMin>0.5){
    //   Double_t ThPhR[3]={DummyRecoCor[0]*(Interferometer::pi/180.),DummyRecoCor[1]*(Interferometer::pi/180.),DummyRecoCor[2]};
    //   Double_t XYZ[3]={0,0,0};
    //   Interferometer::ThPhRtoXYZ(ThPhR,XYZ);
    
    //   //Interferometer::RootThPhR(XYZ,DummyRecoCor,DummyRecoCor,ExpectedUncertainty,ChHitTime,IgnoreCh,ChSNR,DummyMin,Iterations);
    //   Interferometer::MinimizerThPhR(XYZ,DummyRecoCor,DummyRecoCor,ExpectedUncertainty,ChHitTime,IgnoreCh,ChSNR,ChHitOrder,DummyMin,Iterations);
    //   //}
    
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
