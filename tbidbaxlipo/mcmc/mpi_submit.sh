for i in `seq 10`; do
    bsub -a openmpi -W 24:00 -n 13 -q parallel mpirun.lsf /home/jab69/virtualenvs/pysb/bin/python /home/jab69/projects/tBidBaxLipo/tbidbaxlipo/mcmc/nbd_mcmc_mpi_run.py nsteps=1000000 nbd_observables=iBax nbd_sites=c62 random_seed=$i model=taid
done
