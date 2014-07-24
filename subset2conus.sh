#!/bin/bash

if ! [[ -d subset ]]; then mkdir subset; fi

for f in *.nc.gz; do
	echo $f
	gunzip $f
done	

for f in *.nc; do
	echo $f
	ncea -d lat,15.,60. -d lon,200.,300. $f subset/$f
done

