
void PlotRFTxPos(){

   ////Open the file
  std::ifstream ain("RFTxPositions.txt");
  
  int n1=0;////variable for counting total number of data points
  std::string line;
  double dummya[4]={0,0,0,0};////temporary variable for storing data values from the file

  TH3D *h3=new TH3D("","",100,-500,500,100,-500,500,100,-1000,-1);
  
  //Check if file is open and store data
  if(ain.is_open()){
    while (true){
      ain>>dummya[0]>>dummya[1]>>dummya[2]>>dummya[3];
      if( ain.eof() ) break;     
      cout<<n1<<" "<<dummya[0]<<" "<<dummya[1]<<" "<<dummya[2]<<" "<<dummya[3]<<endl;
      
      h3->Fill(dummya[1],dummya[2],-dummya[3]);
      n1++;
    }
  }
  h3->SetMarkerStyle(20);
  h3->Draw();
  
}
