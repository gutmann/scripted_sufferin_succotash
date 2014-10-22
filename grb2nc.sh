#!/usr/bin/env bash
for file in *; do
	echo $file
	ncl_convert2nc $file -e grb -o ncfiles
done
	
	