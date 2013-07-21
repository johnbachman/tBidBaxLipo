for i in `seq 10`; do
    bsub -a openmpi -W 24:00 -n 9 -q parallel mpirun.lsf /home/jab69/virtualenvs/pysb/bin/python /home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/nbd_plate_mcmc_mpi_run.py nsteps=500000 nbd_site=c175 num_confs=3 replicate=2 random_seed=$i
    bsub -a openmpi -W 24:00 -n 9 -q parallel mpirun.lsf /home/jab69/virtualenvs/pysb/bin/python /home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/nbd_plate_mcmc_mpi_run.py nsteps=500000 nbd_site=c179 num_confs=3 replicate=2 random_seed=$i
    bsub -a openmpi -W 24:00 -n 9 -q parallel mpirun.lsf /home/jab69/virtualenvs/pysb/bin/python /home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/nbd_plate_mcmc_mpi_run.py nsteps=500000 nbd_site=c5 num_confs=3 replicate=2 random_seed=$i
    bsub -a openmpi -W 24:00 -n 9 -q parallel mpirun.lsf /home/jab69/virtualenvs/pysb/bin/python /home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/nbd_plate_mcmc_mpi_run.py nsteps=500000 nbd_site=c15 num_confs=3 replicate=2 random_seed=$i
done
