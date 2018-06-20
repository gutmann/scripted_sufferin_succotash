#!/bin/bash

exe=`mktemp test.XXXXX.exe`

if [[ $FC == gfortran ]]; then
    flags=-Wall
fi

if [[ $FC == ifort ]]; then
    flags='-warn all'
fi


echo ${FC} $flags $1 -o ${exe}
${FC} $flags $1 -o ${exe}

./${exe}

rm ${exe}
