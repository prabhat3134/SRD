#!/bin/bash

i=202
i1=26
x=1
while [ $x -le 2 ]
do
	mv RP_$i RP_$i1 
	x=$(( $x +1))
	i=$(( $i +1))
	i1=$(( $i1 +1))
done
