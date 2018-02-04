#!/bin/bash

# parameter_set 1: density ratio=3
i1=125
j=9
x=1
while [ $x -le 8 ]
do
	i2=$(( $i1 + $j ))
	nohup ./run_script.sh $i1 $i2  > runscript_$i1_$i2.out &
	i1=$(( $i2 + 1))
	x=$(( $x + 1))
done
