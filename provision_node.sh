#!/bin/bash

# Exit immediately if a command exits with error
set -e

if [ -e $SGE_CLUSTER_NAME ]; then
    echo "Run this only on a StarCluster node."
    exit -1
fi

if [[ $(hostname) == master ]]; then
    # Put any temp/downloaded files into /data
    cd /data
    # Download BioNetGen if we need to
    if [ ! -e /data/BioNetGen-2.2.5-stable ]; then
        curl -O -J http://mmbios.org/index.php/bionetgen-2-2-5-stable/bionetgen-2-2-5-stable-zip?format=raw
        unzip bionetgen-2.2.5-stable.zip
    fi
    # Get the code if we need to
    if [ ! -e /data/tBidBaxLipo ]; then
        git clone https://github.com/johnbachman/tBidBaxLipo.git
    fi
    # Make sure the code is up-to-date
    cd /data/tBidBaxLipo
    git pull
    git checkout master

    # Get the latest emcee from my remote
    if [ ! -e /data/emcee]; then
        git clone https://github.com/johnbachman/emcee.git
    fi
    # Make sure the code is up-to-date
    cd /data/emcee
    git pull
    git checkout master

    # Make the data dir writeable
    chown -R sgeadmin /data
fi

cd /tmp

# Recompile mpi4py against openmpi
pip uninstall -y mpi4py
update-alternatives --set mpi /usr/lib/openmpi/include
pip install mpi4py

# Install additional software
pip install triangle_plot
pip install ipdb

# Install packages for PySB
pip install sympy
pip install git+https://github.com/pysb/pysb.git
# Put BNG into the right location for this node
cp -a /data/BioNetGen-2.2.5-stable/ /usr/local/share/BioNetGen

# Test PySB
python -m pysb.examples.run_tutorial_a

# Append the following to /etc/apt/sources.list
echo 'deb http://old-releases.ubuntu.com/ubuntu/ raring main' >> /etc/apt/sources.list
echo 'deb-src http://old-releases.ubuntu.com/ubuntu/ raring main' >> /etc/apt/sources.list
echo 'deb http://old-releases.ubuntu.com/ubuntu/ raring-updates main' >> /etc/apt/sources.list
echo 'deb-src http://old-releases.ubuntu.com/ubuntu/ raring-updates main' >> /etc/apt/sources.list
echo 'deb http://old-releases.ubuntu.com/ubuntu/ raring universe' >> /etc/apt/sources.list
echo 'deb-src http://old-releases.ubuntu.com/ubuntu/ raring universe' >> /etc/apt/sources.list
echo 'deb http://old-releases.ubuntu.com/ubuntu/ raring-updates universe' >> /etc/apt/sources.list
echo 'deb-src http://old-releases.ubuntu.com/ubuntu/ raring-updates universe' >> /etc/apt/sources.list

# Update apt
apt-get update
# Install ack
apt-get install ack-grep
# Install texlive (ghostscript installed with texlive)
apt-get install -y texlive
apt-get install -y texlive-latex-extra
apt-get install -y dvipng

