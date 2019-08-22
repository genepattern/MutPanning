#!/bin/bash

MEM=$1
shift
python /mutpanning/adjustmemory.py $MEM > .mem.txt
MEM=$(cat ./.mem.txt)
echo Memory is $MEM


java -Xmx${MEM} -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning $@
JOBDIR=$PWD
cd /mutpanning_temp/SignificanceFiltered/

tar -czvf $JOBDIR/Results.tar.gz *


