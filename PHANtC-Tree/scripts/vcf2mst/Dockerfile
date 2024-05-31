FROM continuumio/miniconda3:4.9.2

RUN conda install -c conda-forge -c bioconda grapetree==2.1 &&\
    apt-get install -y perl


COPY *.pl /usr/local/bin/ 

RUN chmod +x /usr/local/bin/*.pl 

WORKDIR /tmp