import sys
import subprocess

filelist = sys.argv[1:]
for filename in filelist:
    outfile = '%s.evi.out' % filename.split('.')[0]
    errfile = '%s.evi.err' % filename.split('.')[0]
    cmd = 'qsub -b y -cwd -V -o %s -e %s ' \
          'python -m tbidbaxlipo.pt.tabulate_bf %s' % \
          (outfile, errfile, filename)
    cmd_list = cmd.split(' ')
    p = subprocess.Popen(cmd_list)
    out, err = p.communicate()
