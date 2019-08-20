#!/bin/bash
#MAX_SAMPLE=$1
#shift

#SAMPLE_COUNT=$(tr -d -c '\n\r' < $3 | wc -c)
#if  (( $SAMPLE_COUNT > MAX_SAMPLE  )); then
#   >&2 echo "Too many samples, $SAMPLE_COUNT, to process on this GenePattern server."
#   >&2 echo "Because MutPanning is very resource intensive, to process more than $MAX_SAMPLE samples you need to either"
#   >&2 echo "run your own GenePattern server or use the docker container directly to do your analysis"
#   exit 999
#else
#   echo $SAMPLE_COUNT samples is below the max limit of $MAX_SAMPLE
#fi

MEM=$1

shift
python /mutpanning/adjustmemory.py $MEM > .mem.txt
MEM=$(cat ./.mem.txt)
echo Memory is $MEM


java -Xmx${MEM} -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning $@

tar -czvf Results.tar.gz /mutpanning_temp/SignificanceFiltered/


