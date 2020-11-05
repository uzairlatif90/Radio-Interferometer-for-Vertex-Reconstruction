#include "Interferometer.cc"

void RunInterferometerMultTx(){
  DeclareAntennaConfig();
  
  double ChHitTime[2][TotalAntennasRx]; ////Channel Hit Time
  int IgnoreCh[2][TotalAntennasRx]; ////Channel Hit Time
  double ChDRTime[TotalAntennasRx]; ////Channel Hit Time

  //double XYZInitialValues[9][3]={{1,1,1},{0,0,1},{-1,-1,1},{0,1,1},{1,0,1},{0,-1,1},{-1,0,1},{1,-1,1},{-1,1,1}};
  const int InitialValueNum=3;
  double XYZInitialValues[InitialValueNum][3]={{1,1,1},{0,0,1},{-1,-1,1}};
  double InitialValues[3];
  double RecoXYZValues[9][3];
  double RecofMinValues[9];

  double FinalMinValue=0;
  int FinalMinValueBin=0;
  int countTotal=0;
  int countRecoSuccess=0;
  int countRecoFailed=0;
  
  double FinalTxCor[3];
  double JitterNumber=5;

  TH1D *hx=new TH1D("","",200,-500,500);
  TH1D *hy=new TH1D("","",200,-500,500);
  TH1D *hz=new TH1D("","",200,-500,500);

  TH1D *hxZoom=new TH1D("","",100,-50,50);
  TH1D *hyZoom=new TH1D("","",100,-50,50);
  TH1D *hzZoom=new TH1D("","",100,-50,50);
  
  ofstream aout("RFTxPositions.txt");
  
  for(double i=-500; i<501;i=i+100){
    for(double j=-500; j<501;j=j+100){

      if((i!=0 && j==0) || (i==0 && j!=0) || (i!=0 && j!=0)){
	
	for(double k=-1; k>-1001;k=k-50){
	  double DummyTx[3]={i,j,-k};
	  for(int iRx=0;iRx<TotalAntennasRx;iRx++){
	    for(int iray=0;iray<3;iray++){
	      IgnoreCh[iray][iRx]=1;
	    }
	  }
	
	  Interferometer::GenerateChHitTimeAndCheckHits(DummyTx,ChHitTime,IgnoreCh);
	  bool CheckStationTrigger=Interferometer::CheckTrigger(IgnoreCh);
	  if(CheckStationTrigger==true){
	    int status;
	    //Interferometer::ReadChHitTimeFromData("ChHitTimesFromData.txt",ChHitTime);
	    //Interferometer::AddGaussianJitterToHitTimes(JitterNumber,ChHitTime);
	    Interferometer::FindFirstHitAndNormalizeHitTime(ChHitTime,IgnoreCh,ChDRTime);
	   
	    int iInVal=0;
	    while(iInVal<InitialValueNum){
	      InitialValues[0]=XYZInitialValues[iInVal][0];
	      InitialValues[1]=XYZInitialValues[iInVal][1];
	      InitialValues[2]=XYZInitialValues[iInVal][2];
	      Interferometer::Minimizer(InitialValues,RecoXYZValues[iInVal],ChHitTime,IgnoreCh,RecofMinValues[iInVal]);
	      //cout<<"Attempt No. "<<iInVal<<" :: Reco Results are: Xreco="<<RecoXYZValues[iInVal][0]<<" ,Yreco="<<RecoXYZValues[iInVal][1]<<" ,Zreco="<<RecoXYZValues[iInVal][2]<<" "<<RecofMinValues[iInVal]<<endl;
	      if(iInVal==0){
	      	FinalMinValue=RecofMinValues[iInVal];
		FinalMinValueBin=iInVal;
	      }
	      if(RecofMinValues[iInVal]< RecofMinValues[iInVal-1] && RecofMinValues[iInVal]<RecofMinValues[0] && iInVal>0){
		FinalMinValue=RecofMinValues[iInVal];
		FinalMinValueBin=iInVal;
	      }
	      iInVal++;
	    }
	    
	    FinalTxCor[0]=RecoXYZValues[FinalMinValueBin][0];
	    FinalTxCor[1]=RecoXYZValues[FinalMinValueBin][1];
	    FinalTxCor[2]=RecoXYZValues[FinalMinValueBin][2];
	    
	    double CheckReco=0;
	    CheckReco=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
	    cout<<"Final Reco Results are: dX="<<i-FinalTxCor[0]<<" ,dY="<<j-FinalTxCor[1]<<" ,dZ="<<k-FinalTxCor[2]<<" |  Xtrue="<<i<<" ,Ytrue="<<j<<" ,Ztrue="<<k<<" | Xreco="<<FinalTxCor[0]<<" ,Yreco="<<FinalTxCor[1]<<" ,Zreco="<<FinalTxCor[2]<<" "<<FinalMinValue<<endl;
	    if(CheckReco>=3){
	      cout<<"Minimization Failed!"<<endl;
	      aout<<countRecoFailed<<" "<<DummyTx[0]<<" "<<DummyTx[1]<<" "<<DummyTx[2]<<endl;
	      countRecoFailed++;
	    }
	    hx->Fill(i-FinalTxCor[0]);	
	    hy->Fill(j-FinalTxCor[1]);
	    hz->Fill(k-FinalTxCor[2]);

	    hxZoom->Fill(i-FinalTxCor[0]);	
	    hyZoom->Fill(j-FinalTxCor[1]);
	    hzZoom->Fill(k-FinalTxCor[2]);
	    
	    CheckReco=fabs(i-FinalTxCor[0])+fabs(j-FinalTxCor[1])+fabs(k-FinalTxCor[2]);
	    if(CheckReco<3){
	      countRecoSuccess++;
	    }
	 
	    countTotal++;
	  }////Check Station Trigger
	}////k loop
      }////check Tx is not at zero
    }////j loop
  }////i loop

  cout<<"SuccessReco: "<<countRecoSuccess<<", FailedReco: "<<countRecoFailed<<", TotalReco: "<<countTotal<<endl;

  hx->SetTitle("#DeltaX; X_{true} - X_{reco} (m);No. of Events");
  hy->SetTitle("#DeltaY; Y_{true} - Y_{reco} (m);No. of Events");
  hz->SetTitle("#DeltaZ; Z_{true} - z_{reco} (m);No. of Events");
  TCanvas *c1=new TCanvas("c1","c1");
  c1->Divide(3);
  c1->cd(1);
  hx->Draw();
  c1->cd(2);
  hy->Draw();
  c1->cd(3);
  hz->Draw();

  hxZoom->SetTitle("#DeltaX; X_{true} - X_{reco} (m);No. of Events");
  hyZoom->SetTitle("#DeltaY; Y_{true} - Y_{reco} (m);No. of Events");
  hzZoom->SetTitle("#DeltaZ; Z_{true} - z_{reco} (m);No. of Events");
  TCanvas *c12=new TCanvas("c12","c12");
  c12->Divide(3);
  c12->cd(1);
  hxZoom->Draw();
  c12->cd(2);
  hyZoom->Draw();
  c12->cd(3);
  hzZoom->Draw();
 
  
}
