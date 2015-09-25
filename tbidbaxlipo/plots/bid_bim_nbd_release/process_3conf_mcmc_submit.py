import sys
import subprocess

filelist = sys.argv[1:]
for filename in filelist:
    outfile = '%s.process_3conf.out' % filename.split('.')[0]
    errfile = '%s.process_3conf.err' % filename.split('.')[0]
    cmd = 'qsub -b y -cwd -V -o %s -e %s python -m ' \
          'tbidbaxlipo.plots.bid_bim_nbd_release.process_3conf_mcmc %s' % \
          (outfile, errfile, filename)
    cmd_list = cmd.split(' ')
    p = subprocess.Popen(cmd_list)
    out, err = p.communicate()

