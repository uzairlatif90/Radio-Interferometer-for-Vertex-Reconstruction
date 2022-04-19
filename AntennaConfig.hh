#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include <chrono>

////Define your antennas
int TotalAntennasRx=0;
vector <double> AntennaCoordRx[3];
double AvgAntennaCoordRx[3]={0,0,0};

void DeclareAntennaConfig(){
  TotalAntennasRx=16;
  
  int NumAntAtOneDepth=4;
  const int NumRings=1;
  double RingRadius[NumRings]={sqrt(10*10+10*10)};
  const int NumDepths=4;
  double AntDepths[NumDepths]={-180,-185,-195,-200};

  int TotalAntennasRx=0;
  for(int idepth=0;idepth<4;idepth++){
    for(int iring=0;iring<NumRings;iring++){
      for(int ipnt=0;ipnt<NumAntAtOneDepth;ipnt++){  
	
	double x,y,z;
	x=(RingRadius[iring]/sqrt(2))*(iring+1);
	y=(RingRadius[iring]/sqrt(2))*(iring+1);
	z=AntDepths[idepth];
	TVector3 XYvector(x,y,0);
	int isample=0;
	
	XYvector.RotateZ(((360.0/NumAntAtOneDepth))*(TMath::Pi()/180.0)*ipnt);
	x=XYvector.X();
	y=XYvector.Y();
	
	if(fabs(x)<0.0001){
	  x=0;
	}
	if(fabs(y)<0.0001){
	  y=0;
	}

	// AntennaCoordRx[TotalAntennasRx][0]=x;
	// AntennaCoordRx[TotalAntennasRx][1]=y;
	// AntennaCoordRx[TotalAntennasRx][2]=z;

	AntennaCoordRx[0].push_back(x);
	AntennaCoordRx[1].push_back(y);
	AntennaCoordRx[2].push_back(z);

	AvgAntennaCoordRx[0]+=x;
	AvgAntennaCoordRx[1]+=y;
	AvgAntennaCoordRx[2]+=z;
	
	TotalAntennasRx++;
	//cout<<idepth<<" "<<iring<<" "<<ipnt<<" "<<x<<" "<<y<<" "<<z<<endl;
      }
    }
  }

  AvgAntennaCoordRx[0]=AvgAntennaCoordRx[0]/TotalAntennasRx;
  AvgAntennaCoordRx[1]=AvgAntennaCoordRx[1]/TotalAntennasRx;
  AvgAntennaCoordRx[2]=AvgAntennaCoordRx[2]/TotalAntennasRx;

    /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    
    AntennaCoordRx[0][ich]-=AvgAntennaCoordRx[0];
    AntennaCoordRx[1][ich]-=AvgAntennaCoordRx[1];
    AntennaCoordRx[2][ich]-=AvgAntennaCoordRx[2];
    
    cout<<ich<<" "<<AntennaCoordRx[0][ich]<<" "<<AntennaCoordRx[1][ich]<<" "<<AntennaCoordRx[2][ich]<<endl;
    //delete [] antLocRx;
  }

  
}

void DeclareAntennaConfigARA(int StationId){
  TotalAntennasRx=16;
  
  AraStationInfo *stationInfo=new AraStationInfo(StationId,2019);
  //AraStationInfo *stationInfo=new AraStationInfo(StationId);
  AraGeomTool *geom = AraGeomTool::Instance();
  
  /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    double *antLocRx=stationInfo->getAntennaInfo(ich)->getLocationXYZ();
    
    AntennaCoordRx[0].push_back(antLocRx[0]);
    AntennaCoordRx[1].push_back(antLocRx[1]);
    AntennaCoordRx[2].push_back(antLocRx[2]);

    AvgAntennaCoordRx[0]+=antLocRx[0];
    AvgAntennaCoordRx[1]+=antLocRx[1];
    AvgAntennaCoordRx[2]+=antLocRx[2];
    
    //cout<<ich<<" "<<AntennaCoordRx[ich][0]<<" "<<AntennaCoordRx[ich][1]<<" "<<AntennaCoordRx[ich][2]<<endl;
    //delete [] antLocRx;
  }
  // delete geom;
  // delete stationInfo;

  AvgAntennaCoordRx[0]=AvgAntennaCoordRx[0]/TotalAntennasRx;
  AvgAntennaCoordRx[1]=AvgAntennaCoordRx[1]/TotalAntennasRx;
  AvgAntennaCoordRx[2]=AvgAntennaCoordRx[2]/TotalAntennasRx;


  /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    
    AntennaCoordRx[0][ich]-=AvgAntennaCoordRx[0];
    AntennaCoordRx[1][ich]-=AvgAntennaCoordRx[1];
    AntennaCoordRx[2][ich]-=AvgAntennaCoordRx[2];
    
    cout<<ich<<" "<<AntennaCoordRx[0][ich]<<" "<<AntennaCoordRx[1][ich]<<" "<<AntennaCoordRx[2][ich]<<endl;
    //delete [] antLocRx;
  }
}

void DeclareAntennaConfigAraSim(){
  TotalAntennasRx=16;
  
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    AntennaCoordRx[0].push_back(0);
    AntennaCoordRx[1].push_back(0);
    AntennaCoordRx[2].push_back(0);
  }
  
  AntennaCoordRx[0][5]=3.58244;
  AntennaCoordRx[1][5]=-9.72894;
  AntennaCoordRx[2][5]=-189.4;
  AntennaCoordRx[0][13]=3.5876;
  AntennaCoordRx[1][13]=-9.72378;
  AntennaCoordRx[2][13]=-186.115;
  AntennaCoordRx[0][1]=3.61241;
  AntennaCoordRx[1][1]=-9.69901;
  AntennaCoordRx[2][1]=-170.347;
  AntennaCoordRx[0][9]=3.617;
  AntennaCoordRx[1][9]=-9.69443;
  AntennaCoordRx[2][9]=-167.428;
  AntennaCoordRx[0][6]=-3.8531;
  AntennaCoordRx[1][6]=10.0435;
  AntennaCoordRx[2][6]=-191.242;
  AntennaCoordRx[0][14]=-3.84725;
  AntennaCoordRx[1][14]=10.0494;
  AntennaCoordRx[2][14]=-187.522;
  AntennaCoordRx[0][2]=-3.8222;
  AntennaCoordRx[1][2]=10.0745;
  AntennaCoordRx[2][2]=-171.589;
  AntennaCoordRx[0][10]=-3.8173;
  AntennaCoordRx[1][10]=10.0794;
  AntennaCoordRx[2][10]=-168.468;
  AntennaCoordRx[0][7]=-9.11733;
  AntennaCoordRx[1][7]=-3.39657;
  AntennaCoordRx[2][7]=-194.266;
  AntennaCoordRx[0][15]=-9.11217;
  AntennaCoordRx[1][15]=-3.39141;
  AntennaCoordRx[2][15]=-190.981;
  AntennaCoordRx[0][3]=-9.08766;
  AntennaCoordRx[1][3]=-3.36688;
  AntennaCoordRx[2][3]=-175.377;
  AntennaCoordRx[0][11]=-9.08301;
  AntennaCoordRx[1][11]=-3.36223;
  AntennaCoordRx[2][11]=-172.42;
  AntennaCoordRx[0][4]=9.31774;
  AntennaCoordRx[1][4]=3.01173;
  AntennaCoordRx[2][4]=-189.502;
  AntennaCoordRx[0][12]=9.32239;
  AntennaCoordRx[1][12]=3.01638;
  AntennaCoordRx[2][12]=-186.546;
  AntennaCoordRx[0][0]=9.34805;
  AntennaCoordRx[1][0]=3.04201;
  AntennaCoordRx[2][0]=-170.247;
  AntennaCoordRx[0][8]=9.35239;
  AntennaCoordRx[1][8]=3.04635;
  AntennaCoordRx[2][8]=-167.492;
  AvgAntennaCoordRx[0]=0;
  AvgAntennaCoordRx[1]=0;
  AvgAntennaCoordRx[2]=-179.934;

  /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    
    AntennaCoordRx[0][ich]-=AvgAntennaCoordRx[0];
    AntennaCoordRx[1][ich]-=AvgAntennaCoordRx[1];
    AntennaCoordRx[2][ich]-=AvgAntennaCoordRx[2];
    
    cout<<ich<<" "<<AntennaCoordRx[0][ich]<<" "<<AntennaCoordRx[1][ich]<<" "<<AntennaCoordRx[2][ich]<<endl;
    //delete [] antLocRx;
  }

}

void DeclareAntennaConfigRNOG(){
  TotalAntennasRx=24;
  
  AntennaCoordRx[0].resize(TotalAntennasRx);
  AntennaCoordRx[1].resize(TotalAntennasRx);
  AntennaCoordRx[2].resize(TotalAntennasRx);
  
  //// "phased array",
  AntennaCoordRx[0][0]= 0.0;
  AntennaCoordRx[1][0]= 20.0;
  AntennaCoordRx[2][0]= -97.0;

  //// "phased array";
  AntennaCoordRx[0][1]= 0.0;
  AntennaCoordRx[1][1]= 20.0;
  AntennaCoordRx[2][1]= -96.0;
    
  //// "phased array";
  AntennaCoordRx[0][2]= 0.0;
  AntennaCoordRx[1][2]= 20.0;
  AntennaCoordRx[2][2]= -95.0;
    
  //// "phased array";
  AntennaCoordRx[0][3]= 0.0;
  AntennaCoordRx[1][3]= 20.0;
  AntennaCoordRx[2][3]= -94.0;
    
  //// "power string Hpol";
  AntennaCoordRx[0][4]= 0.0;
  AntennaCoordRx[1][4]= 20.0;
  AntennaCoordRx[2][4]= -93.0;
    
  //// "power string Hpol";
  AntennaCoordRx[0][5]= 0.0;
  AntennaCoordRx[1][5]= 20.0;
  AntennaCoordRx[2][5]= -92.0;
    
  //// "power string Vpol";
  AntennaCoordRx[0][6]= 0.0;
  AntennaCoordRx[1][6]= 20.0;
  AntennaCoordRx[2][6]= -80.0;
    
  //// "power string Vpol";
  AntennaCoordRx[0][7]= 0.0;
  AntennaCoordRx[1][7]= 20.0;
  AntennaCoordRx[2][7]= -60.0;
    
  //// "power string Vpol";
  AntennaCoordRx[0][8]= 0.0;
  AntennaCoordRx[1][8]= 20.0;
  AntennaCoordRx[2][8]= -40.0;
    
  //// "Helper string B";
  AntennaCoordRx[0][9]= -17.3;
  AntennaCoordRx[1][9]= -10.0;
  AntennaCoordRx[2][9]= -96.0;
    
  //// "Helper string B";
  AntennaCoordRx[0][10]= -17.3;
  AntennaCoordRx[1][10]= -10.0;
  AntennaCoordRx[2][10]= -95.0;
    
  //// "Helper string B";
  AntennaCoordRx[0][11]= -17.3;
  AntennaCoordRx[1][11]= -10.0;
  AntennaCoordRx[2][11]= -94.0;
    
  //// "first arm north from DAQ; west postition in arm; pointing side-downwards";
  AntennaCoordRx[0][12]= 1.5;
  AntennaCoordRx[1][12]= -11;
  AntennaCoordRx[2][12]= -2;
    
  //// "first arm north from DAQ; middle position in arm; pointing upwards";
  AntennaCoordRx[0][13]= 0.0;
  AntennaCoordRx[1][13]= -11;
  AntennaCoordRx[2][13]= -2.0;
    
  //// "first arm north from DAQ; east position in arm; pointing side-downwards";
  AntennaCoordRx[0][14]= -1.5;
  AntennaCoordRx[1][14]= -11;
  AntennaCoordRx[2][14]= -2;
    
  //// "second arm south-east from DAQ; north-east position in arm; pointing side-downwards";
  AntennaCoordRx[0][15]= -10.276;
  AntennaCoordRx[1][15]= 4.2;
  AntennaCoordRx[2][15]= -2.0;
    
  //// "second arm south-east from DAQ; middle position in arm; pointing upwards";
  AntennaCoordRx[0][16]= -9.526;
  AntennaCoordRx[1][16]= 5.5;
  AntennaCoordRx[2][16]= -2.0;
    
  //// "second arm south-east from DAQ; south-west position in arm; pointing side-downwards";
  AntennaCoordRx[0][17]= -8.776;
  AntennaCoordRx[1][17]= 6.8;
  AntennaCoordRx[2][17]= -2;
    
  //// "third arm south-west from DAQ; south-east position in arm; pointing side-downwards";
  AntennaCoordRx[0][18]= 8.776;
  AntennaCoordRx[1][18]= 6.8;
  AntennaCoordRx[2][18]= -2.0;
    
  //// "third arm south-west from DAQ; middle position in arm; pointing upwards ";
  AntennaCoordRx[0][19]= 9.526;
  AntennaCoordRx[1][19]= 5.5;
  AntennaCoordRx[2][19]= -2.0;
    
  //// "third arm south-west from DAQ; north-west position in arm; pointing side-downwards";
  AntennaCoordRx[0][20]= 10.276;
  AntennaCoordRx[1][20]= 4.2;
  AntennaCoordRx[2][20]= -2.0;
    
  //// "Helper string A";
  AntennaCoordRx[0][21]= 17.3;
  AntennaCoordRx[1][21]= -10.0;
  AntennaCoordRx[2][21]= -96.0;
    
  //// "Helper string B";
  AntennaCoordRx[0][22]= 17.3;
  AntennaCoordRx[1][22]= -10.0;
  AntennaCoordRx[2][22]= -95.0;
    
  //// "Helper string B";
  AntennaCoordRx[0][23]= 17.3;
  AntennaCoordRx[1][23]= -10.0;
  AntennaCoordRx[2][23]= -94.0;
  
  /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    AvgAntennaCoordRx[0]+=AntennaCoordRx[0][ich];
    AvgAntennaCoordRx[1]+=AntennaCoordRx[1][ich];
    AvgAntennaCoordRx[2]+=AntennaCoordRx[2][ich];    
  }

  AvgAntennaCoordRx[0]=AvgAntennaCoordRx[0]/TotalAntennasRx;
  AvgAntennaCoordRx[1]=AvgAntennaCoordRx[1]/TotalAntennasRx;
  AvgAntennaCoordRx[2]=AvgAntennaCoordRx[2]/TotalAntennasRx;

  /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    AntennaCoordRx[0][ich]-=AvgAntennaCoordRx[0];
    AntennaCoordRx[1][ich]-=AvgAntennaCoordRx[1];
    AntennaCoordRx[2][ich]-=AvgAntennaCoordRx[2];
    
    cout<<ich<<" "<<AntennaCoordRx[0][ich]<<" "<<AntennaCoordRx[1][ich]<<" "<<AntennaCoordRx[2][ich]<<endl;
  }
    
}

