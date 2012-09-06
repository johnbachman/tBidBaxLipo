from nbd_data import *

csv_file = open('nbd_data.csv', 'w')
for i, timept in enumerate(time_other):
    csv_file.write('%d, %f, %f, %f, %f, %f\n' % \
            (timept, nbd3c[1][i], nbd62c[1][i], nbd120c[1][i], nbd122c[1][i], nbd126c[1][i]))
csv_file.close()

