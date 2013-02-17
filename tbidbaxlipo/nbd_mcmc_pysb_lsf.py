import sys
import subprocess

num_chains = 20
queue = 'sorger_1d'

if __name__ == '__main__':
    for i in range(0, num_chains):
        cmd_list = ["bsub", "-q", queue,
                    "python", "nbd_mcmc_pysb.py",
                    "random_seed=%d" % i]
        cmd_list += sys.argv[1:]
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)
