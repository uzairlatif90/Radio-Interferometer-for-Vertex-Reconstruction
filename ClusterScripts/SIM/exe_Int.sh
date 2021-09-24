#!/bin/bash
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh
root -q -b -l '/data/user/ulatif/SIM_Inter/ReconstructSIM.C('$1'',''$2'',''$3'',''$4')'
#root -q -b -l '/data/user/ulatif/templatematch/GetBrianEvents.C("'$1'"',''$2'',''$3')'
