In order to run the code you will have to include AraSim's data directory in this folder too

The main script is called ReconstructAraSimevents.C. The script runs over all the events in a given AraSim file, does reconstruction and stores everthing in an output root file. It can be run through an executable "exe_int.sh" by running the following command :

bash exe_Int.sh 400315 /data/user/ulatif/AraSimReco/E2_nomag_new/AraOut.setup_nue_cc_art_E2.txt.run400315.root ./output/Run400315.root

Here the first input is the run number, the second input is the address of the file AraSim file and the third input is the destination and name of the output file where the results will be stored. The executable and the script are located in the same folder. Everything required by the script is in this one folder. So you can copy paste the whole folder into a folder of your own and run from there.
