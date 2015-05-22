for (( i=3; i <= 8; i++ )); do
    mpirun -n $i -x PYTHONPATH -hostfile hostfile_1node python mpi_speed_test.py
done
