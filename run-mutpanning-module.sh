#!/bin/bash

$@
JAVA_OPTS=$1
shift

java $JAVA_OPTS -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning $@

#  java -Xmx55G -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning /opt/gpbeta/gp_home/jobResults/<job_id> <maf.file>  <sample.annotation.file>  <min.samples.per.cluster>  <min.mutations.per.cluster> <job.cpuCount>  <min.cbase.samples> <min.cbase.mutations>  /mutpanning/src/ComputeDistribution.py /mutpanning/Hg19/



tar -czvf Results.tar.gz SignificanceFiltered/a

# remove intermediate files
rm -rf MutationRateClusters EntityCounts PostSignFilter CBASE AlignHg19 AffinityCounts


