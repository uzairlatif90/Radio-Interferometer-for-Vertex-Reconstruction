void MakeEventFile(){

  ofstream aout("SIM_EVENTS.txt");
  int eventcount=0;
  for(int ith=1;ith<180;ith=ith+5){
    for(int iph=0;iph<360;iph=iph+20){
      for(int ir=1;ir<21;ir++){

	cout<<eventcount<<" "<<ith<<" "<<iph<<" "<<ir*50<<endl;
	aout<<eventcount<<" "<<ith<<" "<<iph<<" "<<ir*50<<endl;
	
	eventcount++;
      }
    }
  }
}
