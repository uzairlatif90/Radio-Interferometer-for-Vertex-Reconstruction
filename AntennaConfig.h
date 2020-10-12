////Define your antennas
const Int_t TotalAntennasRx=8;
Double_t AntennaCoordRx[TotalAntennasRx][3];

void DeclareAntennaConfig(){
  AntennaCoordRx[0][0]=10;
  AntennaCoordRx[0][1]=10;
  AntennaCoordRx[0][2]=-5;
  
  AntennaCoordRx[1][0]=10;
  AntennaCoordRx[1][1]=-10;
  AntennaCoordRx[1][2]=-5;
  
  AntennaCoordRx[2][0]=-10;
  AntennaCoordRx[2][1]=10;
  AntennaCoordRx[2][2]=-5;
  
  AntennaCoordRx[3][0]=-10;
  AntennaCoordRx[3][1]=-10;
  AntennaCoordRx[3][2]=-5;
  
  AntennaCoordRx[4][0]=10;
  AntennaCoordRx[4][1]=10;
  AntennaCoordRx[4][2]=-2;
  
  AntennaCoordRx[5][0]=10;
  AntennaCoordRx[5][1]=-10;
  AntennaCoordRx[5][2]=-2;
  
  AntennaCoordRx[6][0]=-10;
  AntennaCoordRx[6][1]=10;
  AntennaCoordRx[6][2]=-2;
  
  AntennaCoordRx[7][0]=-10;
  AntennaCoordRx[7][1]=-10;
  AntennaCoordRx[7][2]=-2;
}
