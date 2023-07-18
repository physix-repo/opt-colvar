#!/bin/bash

rep=$1

X=$(pwd)

cp optCol2.x colvar_x1 input input.template restart.sh optle cor.x weight-change.sh rep_${rep} 

cd rep_${rep}

${X}/rep_${rep}/optCol2.x &> out2

cd ${X}


