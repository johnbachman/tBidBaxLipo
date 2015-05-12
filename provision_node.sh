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
    curl -O -J http://mmbios.org/index.php/bionetgen-2-2-5-stable/bionetgen-2-2-5-stable-zip?format=raw
    unzip bionetgen-2.2.5-stable.zip
    # Get project code
    git clone https://github.com/johnbachman/tBidBaxLipo.git
    # FIXME FIXME
    cd /data/tBidBaxLipo
    git checkout starcluster
    #pip install -e tBidBaxLipo
    # FIXME FIXME
    chown -R sgeadmin /data
fi

cd /tmp

# Recompile mpi4py against openmpi
pip uninstall -y mpi4py
update-alternatives --set mpi /usr/lib/openmpi/include
pip install mpi4py

# Install additional software
pip install emcee

# Install packages for PySB
pip install sympy
pip install git+https://github.com/pysb/pysb.git
# Put BNG into the right location for this node
cp -a /data/BioNetGen-2.2.5-stable/ /usr/local/share/BioNetGen

# Test PySB
python -m pysb.examples.run_tutorial_a


