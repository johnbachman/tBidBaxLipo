"""A script to dump the labeling ratios into a .csv table for inclusion
in the docs."""

import sys
import csv
from tbidbaxlipo.data.parse_bid_bim_nbd_release import labeling_ratios as lr

# Sort the residues in numerical (not alphanumeric) order
sorted_residues = sorted(lr.keys(), key=lambda x: int(x))

if __name__ == '__main__':

    # Check that we got a filename arg
    if len(sys.argv) < 2:
        print "Usage: %s csv_filename" % sys.argv[0]
        sys.exit(1)
    csv_filename = sys.argv[1]

    with open(csv_filename, 'w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"')
        for res in sorted_residues:
            csv_writer.writerow(['NBD-C%s-Bax' % res, '%.2f' % lr[res]])




