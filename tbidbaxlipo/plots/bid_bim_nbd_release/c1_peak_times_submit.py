import sys
import subprocess

filelist = sys.argv[1:]
for filename in filelist:
    outfile = '%s.c1_peak.out' % filename.split('.')[0]
    errfile = '%s.c1_peak.err' % filename.split('.')[0]
    cmd = 'qsub -b y -cwd -V -o %s -e %s ' \
          'python -m tbidbaxlipo.plots.bid_bim_nbd_release.c1_peak_times %s' % \
          (outfile, errfile, filename)
    cmd_list = cmd.split(' ')
    p = subprocess.Popen(cmd_list)
    out, err = p.communicate()

