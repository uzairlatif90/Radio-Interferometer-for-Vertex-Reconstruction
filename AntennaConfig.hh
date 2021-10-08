#include "TRandom3.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
// #include "AraStationInfo.h"
// #include "AraGeomTool.h"
#include <chrono>

////Define your antennas
const Int_t TotalAntennasRx=16;
Double_t AntennaCoordRx[TotalAntennasRx][3];
Double_t AvgAntennaCoordRx[3]={0,0,0};

void DeclareAntennaConfig(){

  Int_t NumAntAtOneDepth=4;
  const Int_t NumRings=1;
  Double_t RingRadius[NumRings]={sqrt(10*10+10*10)};
  const Int_t NumDepths=4;
  Double_t AntDepths[NumDepths]={-180,-185,-195,-200};

  Int_t TotalAntennas=0;
  for(Int_t idepth=0;idepth<4;idepth++){
    for(Int_t iring=0;iring<NumRings;iring++){
      for(Int_t ipnt=0;ipnt<NumAntAtOneDepth;ipnt++){  
	
	Double_t x,y,z;
	x=(RingRadius[iring]/sqrt(2))*(iring+1);
	y=(RingRadius[iring]/sqrt(2))*(iring+1);
	z=AntDepths[idepth];
	TVector3 XYvector(x,y,0);
	Int_t isample=0;
	
	XYvector.RotateZ(((360.0/NumAntAtOneDepth))*(TMath::Pi()/180.0)*ipnt);
	x=XYvector.X();
	y=XYvector.Y();
	
	if(fabs(x)<0.0001){
	  x=0;
	}
	if(fabs(y)<0.0001){
	  y=0;
	}

	AntennaCoordRx[TotalAntennas][0]=x;
	AntennaCoordRx[TotalAntennas][1]=y;
	AntennaCoordRx[TotalAntennas][2]=z;

	AvgAntennaCoordRx[0]+=x;
	AvgAntennaCoordRx[1]+=y;
	AvgAntennaCoordRx[2]+=z;
	
	TotalAntennas++;
	cout<<idepth<<" "<<iring<<" "<<ipnt<<" "<<x<<" "<<y<<" "<<z<<endl;
      }
    }
  }

  AvgAntennaCoordRx[0]=AvgAntennaCoordRx[0]/TotalAntennasRx;
  AvgAntennaCoordRx[1]=AvgAntennaCoordRx[1]/TotalAntennasRx;
  AvgAntennaCoordRx[2]=AvgAntennaCoordRx[2]/TotalAntennasRx;
  
}

void DeclareAntennaConfigARA(Int_t StationId){
  
  AraStationInfo *stationInfo=new AraStationInfo(StationId,2019);
  //AraStationInfo *stationInfo=new AraStationInfo(StationId);
  AraGeomTool *geom = AraGeomTool::Instance();
  
  /*Get the coordinates of the station channels*/
  for(int ich=0;ich<TotalAntennasRx;ich++) {
    Double_t *antLocRx=stationInfo->getAntennaInfo(ich)->getLocationXYZ();
    
    AntennaCoordRx[ich][0]=antLocRx[0];
    AntennaCoordRx[ich][1]=antLocRx[1];
    AntennaCoordRx[ich][2]=antLocRx[2];

    AvgAntennaCoordRx[0]+=antLocRx[0];
    AvgAntennaCoordRx[1]+=antLocRx[1];
    AvgAntennaCoordRx[2]+=antLocRx[2];
    
    cout<<ich<<" "<<AntennaCoordRx[ich][0]<<" "<<AntennaCoordRx[ich][1]<<" "<<AntennaCoordRx[ich][2]<<endl;
    //delete [] antLocRx;
  }
  // delete geom;
  // delete stationInfo;

  AvgAntennaCoordRx[0]=AvgAntennaCoordRx[0]/TotalAntennasRx;
  AvgAntennaCoordRx[1]=AvgAntennaCoordRx[1]/TotalAntennasRx;
  AvgAntennaCoordRx[2]=AvgAntennaCoordRx[2]/TotalAntennasRx;
}

void DeclareAntennaConfigAraSim(){
  AntennaCoordRx[5][0]=3.58244;
  AntennaCoordRx[5][1]=-9.72894;
  AntennaCoordRx[5][2]=-189.4;
  AntennaCoordRx[13][0]=3.5876;
  AntennaCoordRx[13][1]=-9.72378;
  AntennaCoordRx[13][2]=-186.115;
  AntennaCoordRx[1][0]=3.61241;
  AntennaCoordRx[1][1]=-9.69901;
  AntennaCoordRx[1][2]=-170.347;
  AntennaCoordRx[9][0]=3.617;
  AntennaCoordRx[9][1]=-9.69443;
  AntennaCoordRx[9][2]=-167.428;
  AntennaCoordRx[6][0]=-3.8531;
  AntennaCoordRx[6][1]=10.0435;
  AntennaCoordRx[6][2]=-191.242;
  AntennaCoordRx[14][0]=-3.84725;
  AntennaCoordRx[14][1]=10.0494;
  AntennaCoordRx[14][2]=-187.522;
  AntennaCoordRx[2][0]=-3.8222;
  AntennaCoordRx[2][1]=10.0745;
  AntennaCoordRx[2][2]=-171.589;
  AntennaCoordRx[10][0]=-3.8173;
  AntennaCoordRx[10][1]=10.0794;
  AntennaCoordRx[10][2]=-168.468;
  AntennaCoordRx[7][0]=-9.11733;
  AntennaCoordRx[7][1]=-3.39657;
  AntennaCoordRx[7][2]=-194.266;
  AntennaCoordRx[15][0]=-9.11217;
  AntennaCoordRx[15][1]=-3.39141;
  AntennaCoordRx[15][2]=-190.981;
  AntennaCoordRx[3][0]=-9.08766;
  AntennaCoordRx[3][1]=-3.36688;
  AntennaCoordRx[3][2]=-175.377;
  AntennaCoordRx[11][0]=-9.08301;
  AntennaCoordRx[11][1]=-3.36223;
  AntennaCoordRx[11][2]=-172.42;
  AntennaCoordRx[4][0]=9.31774;
  AntennaCoordRx[4][1]=3.01173;
  AntennaCoordRx[4][2]=-189.502;
  AntennaCoordRx[12][0]=9.32239;
  AntennaCoordRx[12][1]=3.01638;
  AntennaCoordRx[12][2]=-186.546;
  AntennaCoordRx[0][0]=9.34805;
  AntennaCoordRx[0][1]=3.04201;
  AntennaCoordRx[0][2]=-170.247;
  AntennaCoordRx[8][0]=9.35239;
  AntennaCoordRx[8][1]=3.04635;
  AntennaCoordRx[8][2]=-167.492;
  AvgAntennaCoordRx[0]=0;
  AvgAntennaCoordRx[1]=0;
  AvgAntennaCoordRx[2]=-179.934;
}
