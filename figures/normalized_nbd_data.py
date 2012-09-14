"""
Generates a PDF report figure in which each of the NBD trajectories are aligned
on top of each other by fitting their min/max normalization constants to
minimize differences.
"""

from data.nbd_data import nbd120c, nbd122c, nbd126c, nbd3c, time_other as time
from data.nbd_data_62c import nbd62c, time as time_62c
from nbd_analysis import normalize_fit
import matplotlib.pyplot as plt
from numpy import mean

norm_120c = normalize_fit(nbd120c)
norm_122c = normalize_fit(nbd122c)
norm_126c = normalize_fit(nbd126c)
norm_3c = normalize_fit(nbd3c)
norm_62c = normalize_fit(nbd62c)

norm_mean_120c = mean(norm_120c, 0)
norm_mean_122c = mean(norm_122c, 0)
norm_mean_126c = mean(norm_126c, 0)
norm_mean_3c = mean(norm_3c, 0)
norm_mean_62c = mean(norm_62c, 0)

nbd_norm = [norm_120c, norm_122c, norm_126c, norm_3c]
nbd_titles = ['NBD 120c', 'NBD 122c', 'NBD 126c', 'NBD 3c']

plt.ion()

# Plot everything except 62c
for i in range(0, len(nbd_norm)):
    plt.figure()
    plt.title(nbd_titles[i])

    for replicate in nbd_norm[i]:
        plt.plot(time, replicate)

# Plot 62c 
plt.figure()
plt.title('NBD 62c')

for replicate in norm_62c:
    plt.plot(time_62c, replicate)

plt.show()
