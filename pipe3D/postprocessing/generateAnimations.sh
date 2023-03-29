#!/bin/bash

rm -rf animations

mkdir animations

convert -delay 5 rho*.ppm animations/density.gif
convert -delay 5 tau*.ppm animations/shearStress.gif
convert -delay 5 u*.ppm animations/velocity.gif

cd animations

ffmpeg -i density.gif density.avi
ffmpeg -i shearStress.gif shearStress.avi
ffmpeg -i velocity.gif velocity.avi

ffmpeg -i density.avi -i velocity.avi -filter_complex "[0:v:0]pad=(iw*2):ih+1[bg]; [bg][1:v:0]overlay=w" densityVelocity.avi
ffmpeg -i shearStress.avi -i velocity.avi -filter_complex "[0:v:0]pad=(iw*2):ih+1[bg]; [bg][1:v:0]overlay=w" shearStressVelocity.avi
