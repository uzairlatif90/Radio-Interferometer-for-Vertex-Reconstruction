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
