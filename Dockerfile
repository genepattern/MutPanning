FROM ubuntu

RUN mkdir -p /mutpanning

WORKDIR /mutpanning

RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install unzip && \
	apt-get install -y wget openjdk-8-jdk 
	
	
COPY . /mutpanning

RUN javac -cp bin/commons-math3-3.6.1.jar:bin/jdistlib-0.4.5-bin.jar src/AffinityCount_Cosmic.java src/AffinityCount.java src/AlignHG19.java src/CBASE_Solutions.java src/ClusteringEntity.java src/ClusteringPanCancer.java src/ComputeMutationRateClusters_Entities.java src/ComputeSignificance.java src/CountDestructiveMutations.java src/Filter_Step1.java src/Filter_Step2.java src/Filter_Step3.java src/MutPanning.java src/ReformatCBASE.java && \
	mv src/*.class bin

RUN mkdir /mutpanning_temp
RUN mkdir /mutpanning_hg19 && \
    mv /mutpanning/Hg19.zip /mutpanning_hg19 && \
    cd /mutpanning_hg19 && \
    unzip Hg19.zip

RUN apt-get clean && \
    apt-get install -y python 

#CMD ["bash /mutpanning/launch-mutpanning.sh"]

