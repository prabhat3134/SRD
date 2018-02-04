#!/bin/bash

home_dir=$(pwd)
mkdir -p $home_dir/results
cp ../RP_rho $home_dir/
cp ../RP_temp $home_dir/
cp ../RP_vel $home_dir/
nohup ./RP_rho &
nohup ./RP_temp &
nohup ./RP_vel 
cp ../RP_avg $home_dir/
nohup ./RP_avg &
rm -rf RP_* &
