#!/bin/bash
rm -f library_copy.o
rm -f srd_library.mod
home_dir=$(pwd)
dest_dir=$(echo ".")

RP_start=$1
RP_end=$2
x=$RP_start

while [ $x -le $RP_end ] 
do
	# make the directory
	dir_name=$(echo "RP_$x")
	mkdir -p $dest_dir/$dir_name
	
	# copy the executable to the directory
	cp $dest_dir/fish_serial $dest_dir/$dir_name
	cd $dest_dir/$dir_name

	# run the executable 
	nohup ./fish_serial 

	# return to home dir
	cd $home_dir
 
	# update x
	x=$(( $x + 1 )) 
	
	# sleep for generating different set of random numbers
	sleep .5
done
if [ $x -ge 200 ]
then
	sleep 10m
	./avgg.sh
fi
