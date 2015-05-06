# Exit immediately if a command exits with error
set +e

# Put any temp/downloaded files into /home/sgeadmin
cd /home/sgeadmin

# Recompile mpi4py against openmpi
pip uninstall -y mpi4py
update-alternatives --set mpi /usr/lib/openmpi/include
pip install mpi4py

# Install packages for PySB
pip install sympy
pip install git+https://github.com/pysb/pysb.git
curl -O -J http://mmbios.org/index.php/bionetgen-2-2-5-stable/bionetgen-2-2-5-stable-zip?format=raw
unzip bionetgen-2.2.5-stable.zip
mv BioNetGen-2.2.5-stable/ /usr/local/share/BioNetGen

# Test PySB
python -m pysb.examples.run_tutorial_a

# Install additional software
pip install emcee
git clone https://github.com/johnbachman/tBidBaxLipo.git


