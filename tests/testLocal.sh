# quick test of MutPanning in container

# immediately below cached from exec.sh on gp-beta-ami

# NEED TO SET MEMORY FOR XMX TO NUMBER PASSED IN 
# COMMAND_LINE="java -Xmx55G -cp bin/commons-math3-3.6.1.jar:bin/jdistlib-0.4.5-bin.jar:bin MutPanning /mutpanning/"

#  argument 0: root file, where all the other files can be found
#  argument 1: maf file (standard value: root file/MutationsComplete.maf)
#  argument 2: sample annotation file (standard value: root file/SamplesComplete.txt)
#  argument3 : minimal no. samples per cluster (standard value 3)
#  argument4: minimal no. mutations per cluster (standard value 1000)
#  argument5: no cpus available to distribute clustering & CBASE workload (standard value 24)
#  arugment 6: min no. samples for CBASE (standard value 100)
#  argument 7: min no. mutations for CBASE (standard value 5000)
#  argument 8: path to python script to compute distribution of nonsynonymous mutation counts (standard value: root file/src/ComputeDistribution.py)
#  argument 9: path to Hg19 folder (standard value root file/Hg19/)

DATA_DIR=/Users/liefeld/GenePattern/gp_dev/docker/docker-mutpanning/mutpanning_testfiles/Melanoma/
N_CPU=4

COMMAND_LINE="java -Xmx2G -cp /mutpanning/bin/commons-math3-3.6.1.jar:/mutpanning/bin/jdistlib-0.4.5-bin.jar:/mutpanning/bin MutPanning /mutpanning/ $DATA_DIR/MutationsMelanoma.maf $DATA_DIR/SamplesMelanoma.txt 3  1000 $N_CPU 100 5000 /mutpanning/src/ComputeDistribution.py /mutpanning/Hg19/"



LOCAL_DIR=$PWD
ROOT_DIR='/Users/liefeld/GenePattern/gp_dev/docker/docker-mutpanning'

docker run -v $ROOT_DIR:$ROOT_DIR -w $LOCAL_DIR/Job_1 -t genepattern/docker-mutpanning $COMMAND_LINE



