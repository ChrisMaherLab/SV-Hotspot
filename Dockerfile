FROM perl:5.22
#FROM ubuntu:18.04
#FROM r-base

MAINTAINER Abdallah Eteleeb <eteleeb@gmail.com>

LABEL \
    description="SV-HotSpot is a structural variant hotspots detection tool. \
                  It detects SVs and determine their effect on nearby gene expression \
                  using whole-genome sequencing data." \
    version="1.0.0"

##########################################################################################
# Preparation
##########################################################################################
RUN apt-get update && apt-get install -y \
    curl \
    make \ 
    wget \
    cpanminus


RUN ["cpanm", "List::MoreUtils", "Data::Dumper" ] 

##########################################################################################
# INSTALL R
##########################################################################################
RUN apt-get install -y r-base 

##########################################################################################
# INSTALL R PACKAGES
##########################################################################################
RUN R -e "install.packages(c('peakPick', 'ggplot2', 'reshape2', 'gridExtra', 'plyr', 'gtable', 'ggsignif','RCircos', 'data.table'), \
                           repos='http://cran.us.r-project.org')"

##########################################################################################
# INSTALL bedtools v2.28.0
##########################################################################################

WORKDIR /usr/bin
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz 
RUN tar -zxvf bedtools-2.28.0.tar.gz

WORKDIR /usr/bin/bedtools2
RUN make
RUN rm -f /usr/bin/bedtools-2.28.0.tar.gz

WORKDIR /usr/bin/
RUN cp bedtools2/bin/* /usr/bin/

##########################################################################################
###### SV-HotSpot 
##########################################################################################
RUN curl -L https://sourceforge.net/projects/sv-hotspot/files/SV-HotSpot.1.0.0.tar.gz > /usr/bin/SV-HotSpot.tar.gz
RUN tar -xzf /usr/bin/SV-HotSpot.tar.gz -C /usr/bin/

### install sv-hotspot 
WORKDIR /usr/bin
RUN /usr/bin/install.sh -o /usr/bin/
WORKDIR /usr/bin

#remove original files
RUN rm -f /usr/bin/SV-HotSpot.tar.gz

WORKDIR /opt/
COPY Dockerfile /opt/Docerfile

WORKDIR /usr/bin
RUN mv /usr/bin/sv-hotspot.pl /usr/bin/sv-hotspot
RUN mv /usr/bin/plot-peak.pl /usr/bin/plot-peak


