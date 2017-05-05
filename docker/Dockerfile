FROM ubuntu:15.10
MAINTAINER "Greg Way" <gregway@mail.med.upenn.edu>

# Install Ubuntu dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
		git \
        ed \		
		less \
		locales \
		vim-tiny \
		wget \
        ca-certificates \
		bison \
		flex \
		apt-utils \
		gawk \
		libcurl4-openssl-dev \
		libxml2-dev\
		libssl-dev

# Install R
RUN apt-get install -y --no-install-recommends \
		r-base \
		r-base-dev \
		r-recommended \
		r-cran-xml \
		r-cran-rjava

# Install Sleipnir
RUN apt-get install -y --no-install-recommends \
		mercurial \
		gengetopt \
		libboost-regex-dev \
		libboost-graph-dev \
		liblog4cpp5-dev \
		build-essential \
		libgsl0-dev

RUN wget -nv http://download.joachims.org/svm_perf/current/svm_perf.tar.gz \
--directory-prefix svmperf/
RUN cd svmperf && tar -xf svm_perf.tar.gz && make

RUN hg clone https://bitbucket.org/libsleipnir/sleipnir

RUN cd sleipnir && ./gen_auto && ./gen_tools_am
RUN cd sleipnir && ./configure --with-svm-perf=~/svmperf/ && make && make install

# Install hgsc_subtypes
# NOTE - Follow these steps to break cache download following git update
# Because we must allow for the updating of the repository, break this cache
# if the commit ID is different. Run in command line before building image
# export CACHEBUST=`git rev-parse HEAD`
# docker build -t gregway/hgsc_subtypes --build-arg CACHEBUST=$CACHEBUST .
# ARG CACHEBUST=1

# Install git repository
RUN git clone https://github.com/greenelab/hgsc_subtypes.git

# Install curatedOvarianData_1.8.0 from cached download
RUN wget -nv https://bitbucket.org/gwaygenomics/download/raw/\
575d017054fb27a6a0cd62102f86c0ca13bc32a2/ovca/curatedOvarianData_1.8.0.tar.gz
RUN tar -xf curatedOvarianData_1.8.0.tar.gz

