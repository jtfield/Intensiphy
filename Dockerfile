# syntax=docker/dockerfile:1
FROM ubuntu:latest
WORKDIR /project

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=America/New_York

RUN apt-get update && apt-get install -y \ 
    bcftools \
    git \
    tzdata \
    python3 \
    python3-pip \
    raxml \
    seqtk \
    samtools \
    seaview \
    wget \
 && apt-get -y autoremove \
 && apt-get clean autoclean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 
    

#install python and symlink python3 to python
#RUN apt-get install -y python3
#RUN apt-get install -y python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip3 install dendropy --no-cache-dir

#install dependencies available from apt-get
#RUN apt-get install raxml
#RUN apt-get install seqtk
#RUN apt-get install bcftools -y
#RUN apt-get install samtools -y
#RUN apt-get install git -y
#RUN apt-get install seaview -y

#fastx toolkit and bwa-mem need wget to be installed
#RUN apt-get install wget -y
RUN wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 \
 && tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 \
 && cp bin/* /usr/local/bin/ \
 && rm -rf bin fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2  
#RUN tar -xjf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
#RUN cp ./bin/* /usr/local/bin/

RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
 && tar -xjf bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
 && cp ./bwa-mem2-2.0pre2_x64-linux/bwa-mem2* /usr/local/bin/ \
 && rm -rf bwa-mem2-2.0pre2_x64-linux.tar.bz2 bwa-mem2-2.0pre2_x64-linux
#RUN tar -xjf bwa-mem2-2.0pre2_x64-linux.tar.bz2
#RUN cp ./bwa-mem2-2.0pre2_x64-linux/bwa-mem2* /usr/local/bin/

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz
RUN tar -xvzf sratoolkit.3.0.10-ubuntu64.tar.gz
RUN cp -r ./sratoolkit.3.0.10-ubuntu64/bin/* /usr/local/bin
RUN rm -rf sratoolkit.3.0.10-ubuntu64.tar.gz sratoolkit.3.0.10-ubuntu64

# install Extensiphy
RUN git clone https://github.com/McTavishLab/extensiphy.git

# install Intensiphy
RUN git clone https://github.com/jtfield/Intensiphy.git

ENV PATH "$PATH:/project/extensiphy"

WORKDIR /project/Intensiphy

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
