#!/bin/bash

if ! [[ -d subset ]]; then mkdir subset; fi

N=5

let i=0
for f in *.nc.gz; do
	if [[ -e $f ]]; then
		((i=i%N)); ((i++==0)) && wait
		echo $f
		gunzip $f &
	fi
done
wait

subset(){
	ncrcat -d lat,15.,90. -d lon,180.,300. $1 subset/$1
	if [[ -e subset/$1 ]]; then 
		rm $1
	fi
}

let i=0
for f in *.nc; do
	((i=i%N)); ((i++==0)) && wait
	echo $f
	subset $f &
done
