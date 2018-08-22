FROM ubuntu:17.10

RUN mkdir -p /mutpanning

WORKDIR /mutpanning

RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y wget openjdk-8-jdk && \
	apt-get install -y python-pip python2.7 python-dev && \
	pip install scipy==1.1.0 mpmath==1.0.0 
	
COPY . /mutpanning

RUN javac -cp bin/commons-math3-3.6.1.jar:bin/jdistlib-0.4.5-bin.jar src/AffinityCount_Cosmic.java src/AffinityCount.java src/AlignHG19.java src/CBASE_Solutions.java src/ClusteringEntity.java src/ClusteringPanCancer.java src/ComputeMutationRateClusters.java src/ComputeMutationRateEntities.java src/ComputeSignificance_Uniform.java src/ComputeSignificance.java src/CountDestructiveMutations.java src/Filter_Step1.java src/Filter_Step2.java src/Filter_Step3.java src/MutPanning.java src/ReformatCBASE.java && \
	mv src/*.class bin

# 
# Note the Hg19.tar.gz file is not in github but is available at https://datasets.genepattern.org/data/module_support_files/MutPanning/Hg19.tar.gz
#
RUN tar -xvzf Hg19.tar.gz && \
    rm Hg19.tar.gz && \
    chmod -R a+rwx /mutpanning/Hg19

#CMD ["bash /mutpanning/launch-mutpanning.sh"]

