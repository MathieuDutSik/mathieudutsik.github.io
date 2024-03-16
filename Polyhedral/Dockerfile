# ----------------------------------------------------
# Short manual on how to use docker for this test case.
# ----------------------------------------------------
#
# It is used in the following way:
# ---Copy this "Dockerfile" to a directory, e.g. "docker".
#    Do NOT rename the file.
# ---Go to that directory with "cd" command.
# ---Compile the dockerfile with the command
#    docker build --no-cache -t standard_gap_polyhedral:build-env .
#    (the --no-cache command is needed because if the remote source changes
#     this is invisible in the Dockerfile)
# ---Run the docker and access a bash by doing
#    docker run -i -t standard_gap_polyhedral:build-env /bin/bash
#    (This gives you access to a shell where you can run gap for example by "gap.sh")
# ---Examples can be run by e.g. doing
#    cd /opt/gap-4.11.0/pkg/polyhedral/Examples/C-type_Enumeration
#    gap.sh Check_Dim5_Ctype.g


# ----------------------------------------------------
# Compilation of all the commands for polyhedral
# ----------------------------------------------------

# Set the base image
FROM ubuntu:22.04


# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
MAINTAINER Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>


# Install Ubuntu packages
RUN apt-get update -y
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe
RUN apt-get install -y g++ gcc gfortran
RUN apt-get install -y git cmake
RUN apt-get install -y m4
RUN apt-get install -y libboost-all-dev


# GMP
RUN apt-get install -y libgmp-dev
ENV GMP_INCDIR "/usr/include"
ENV GMP_CXX_LINK "-lm -lgmp -lgmpxx"


# GSL
RUN apt-get install -y libgsl-dev


# Required for LINBOX
RUN apt-get install -y autoconf
RUN apt-get install -y libntl-dev
RUN apt-get install -y csh


# EIGEN and /opt directory
WORKDIR /opt
RUN cd /opt && git clone https://gitlab.com/libeigen/eigen.git
ENV EIGEN_PATH /opt/eigen


# mpfrc++
RUN apt-get install -y libmpfr-dev
ENV MPFR_INCDIR /usr/include
ENV MPFR_LINK "-lmpfr"
ENV MPREAL_INCLUDE "-I/opt/mpfrc++-3.6.2"



# Installing gap
RUN wget https://www.gap-system.org/pub/gap/gap-4.11/tar.gz/gap-4.11.0.tar.gz
RUN tar -vxzf gap-4.11.0.tar.gz && rm gap-4.11.0.tar.gz && cd gap-4.11.0 && ./configure && make
RUN cd gap-4.11.0/pkg/grape-4.8.3 && ./configure ../.. && make
RUN cd gap-4.11.0/pkg/nq-2.5.4 && ./configure && make
RUN cd gap-4.11.0/pkg/CaratInterface && ./configure ../.. && make
RUN cd /usr/local/bin && ln -sf /opt/gap-4.11.0/bin/gap.sh


# polyhedral GAP code compilation
RUN mkdir listlink
WORKDIR /opt/listlink


# Downloading the polyhedral code and compiling it
RUN cd /opt/gap-4.11.0/pkg && wget http://mathieudutour.altervista.org/Polyhedral/polyhedral.tar.bz2
RUN cd /opt/gap-4.11.0/pkg && bzip2 -d polyhedral.tar.bz2
RUN cd /opt/gap-4.11.0/pkg && tar -vxf polyhedral.tar && rm polyhedral.tar
RUN cd /opt/gap-4.11.0/pkg/polyhedral && ./configure


# Setting up the gap.ini for the root
RUN mkdir /root/.gap
RUN touch /root/.gap/gap.ini
RUN echo "SetUserPreference( \"InfoPackageLoadingLevel\", 4 );" >> /root/.gap/gap.ini
RUN echo "SetUserPreference( \"PackagesToLoad\", [ \"polyhedral\" ]);" >> /root/.gap/gap.ini
WORKDIR /root
