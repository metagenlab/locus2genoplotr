# Base Image
FROM continuumio/miniconda3:4.8.2
 
################## METADATA ######################
 
LABEL base.image="miniconda3:4.8.2"
LABEL version="1.1"
LABEL software="locus2genoplotr"
LABEL software.version="1.1"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Trestan Pillonel
 
################## INSTALLATION ######################
RUN conda install conda=4.7.12

COPY env.yaml ./
RUN conda env create -f env.yaml
RUN conda clean --all --yes

RUN conda init bash

WORKDIR /data/
ENV PATH /opt/conda/envs/locus2genoplotr/bin:$PATH
