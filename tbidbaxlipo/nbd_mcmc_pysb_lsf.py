import sys
import subprocess

num_chains = 10
queue = 'short'

if __name__ == '__main__':
    # Keyword args are set at the command line as e.g., key=val
    # and subsequently split at the equals sign
    kwargs = dict([arg.split('=') for arg in sys.argv[1:]])

    for i in range(0, num_chains):
        output_filename = '%s_%s_%s_%d_s%d' % \
                (kwargs['model'], kwargs['nbd_site'], kwargs['nbd_observable'],
                 kwargs['nsteps'], i)
        cmd_list = ["bsub", "-q", queue, "-W", "12:00",
                    "-o", output_filename,
                    "python", "nbd_mcmc_pysb.py",
                    "random_seed=%d" % i]
        cmd_list += sys.argv[1:]
        print ' '.join(cmd_list)
        subprocess.call(cmd_list)
