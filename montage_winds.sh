#!/bin/bash

# e.g. montage_winds.sh wfiles/wind_20\*6e5\*.png w20.png
# for i in 5 10 15 20; do  montage_winds.sh wfiles/wind_${i}\*6e5\*.png w${i}.png; done
# montage $1 -crop '1010x554+164+71' -label '%f' -geometry +1+1 -tile 3x4 $2
montage $1 -crop '1050x594+124+61' -label '%f' -geometry +1+1 -tile 3x4 $2