#!/bin/bash
#this line is for loading the correct environment variables in cobalt
source /cvmfs/ara.opensciencegrid.org/trunk/centos7/setup.sh

root -q -b -l 'ReconstructAraSimevents.C('$1'','"'$2'"','"'$3'")'
