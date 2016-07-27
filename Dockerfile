FROM ubuntu
RUN apt-get update
RUN apt-get -f install
RUN apt-get -y install zip unzip
RUN apt-get -y install build-essential
RUN apt-get -y install wget
RUN apt-get install -y python python-pip wget
RUN wget --no-check-certificate https://snap.stanford.edu/releases/Snap-2.4.zip; \
    unzip Snap-2.4.zip; \
    rm -rf Snap-2.4.zip; \
    cd Snap-2.4; make all
