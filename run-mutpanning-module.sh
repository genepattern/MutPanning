#!/bin/bash
MAX_SAMPLE=$1
shift

SAMPLE_COUNT=$(tr -d -c '\n\r' < $3 | wc -c)
if  (( $SAMPLE_COUNT > MAX_SAMPLE  )); then
   >&2 echo "Too many samples, $SAMPLE_COUNT, to process on this GenePattern server."
   >&2 echo "Because MutPanning is very resource intensive, to process more than $MAX_SAMPLE samples you need to either"
   >&2 echo "run your own GenePattern server or use the docker container directly to do your analysis"
   exit 999
else
   echo $SAMPLE_COUNT samples is below the max limit of $MAX_SAMPLE
fi

MEM=$1
shift

echo "RUNNING:   java $JAVA_OPTS -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning $PWD $@"

java $MEM -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning $PWD $@

#  java -Xmx55G -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning /opt/gpbeta/gp_home/jobResults/<job_id> <maf.file>  <sample.annotation.file>  <min.samples.per.cluster>  <min.mutations.per.cluster> <job.cpuCount>  <min.cbase.samples> <min.cbase.mutations>  /mutpanning/src/ComputeDistribution.py /mutpanning/Hg19/



tar -czvf Results.tar.gz SignificanceFiltered/

# remove intermediate files
#tar -czvf IntermediateFiles.tar.gz MutationRateClusters EntityCounts PostSignFilter CBASE AlignHg19 AffinityCounts ClusteringComplete CountDestructive Clustering SignificanceRaw
rm -rf MutationRateClusters EntityCounts PostSignFilter CBASE AlignHg19 AffinityCounts ClusteringComplete CountDestructive Clustering SignificanceFiltered SignificanceRaw


