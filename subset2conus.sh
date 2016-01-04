#!/bin/bash

if ! [[ -d subset ]]; then mkdir subset; fi

for f in *.nc.gz; do
	if [[ -e $f ]]; then
		echo $f
		gunzip $f
	fi
done	

for f in *.nc; do
	echo $f
	ncea -d lat,15.,90. -d lon,180.,300. $f subset/$f
    rm $f
done

