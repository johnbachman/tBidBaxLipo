import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from tbidbaxlipo.util import fitting
import sys

plt.ion()

num_buckets = 100
#num_coins = 1000
max_time = 3
num_reps = 20
p = 0.2
#avg_coins = num_coins / float(num_buckets)

# Analytical integral
def prob_perm(num_coins, num_buckets):
    avg_coins = num_coins / float(num_buckets)
    prob_total = 0
    for k in range(num_coins):
        prob_total += stats.poisson.pmf(k, avg_coins) * \
                      (1 - (1 - p) ** k)
    return prob_total

min_coins = 0
max_coins = 1000
coins_interval = 10

coin_arr = np.arange(min_coins, max_coins + 1, coins_interval)
prob_arr = np.zeros(len(coin_arr))
fmax_arr = np.zeros(len(coin_arr))
k_arr_num = np.zeros(len(coin_arr))
k_arr_anal = np.zeros(len(coin_arr))

def k_func(avg_coins):
    return np.exp(-avg_coins * p) * (-1 + np.exp(avg_coins * p))

for i, coins in enumerate(coin_arr):
    avg_coins = coins / float(num_buckets)
    prob_arr[i] = prob_perm(coins, num_buckets)
    fmax_arr[i] = (1 - np.exp(-avg_coins))
    k_arr_num[i] = prob_arr[i] / float(fmax_arr[i])
    k_arr_anal[i] = k_func(avg_coins) / float(fmax_arr[i])

plt.figure()
plt.plot(coin_arr, prob_arr, color='b', label='Initial prob')
plt.plot(coin_arr, fmax_arr, color='r', label='Poisson Fmax')
plt.plot(coin_arr, k_arr_num, color='g', label='Numerical k')
plt.plot(coin_arr, k_arr_anal, color='m', label='Analytical k')
plt.legend(loc='lower right')
sys.exit()

rep_matrix = np.zeros((num_reps, max_time))

for rep_index in range(num_reps):
    buckets = np.zeros(num_buckets)
    bucket_status = np.ones(num_buckets)
    timecourse = np.zeros(max_time)
    bucket_assignments = np.random.randint(num_buckets, size=num_coins)
    for bucket_index in bucket_assignments:
        buckets[bucket_index] += 1
    num_perm = 0
    for i in range(1, max_time):
        # update buckets
        num_perm = num_buckets - np.sum(bucket_status)
        for bucket_index in np.where(bucket_status == 1)[0]:
            prob_intact = stats.binom.pmf(0, buckets[bucket_index], p)
            # If our random number is greater than prob intact, we permeabilize
            if np.random.rand() > prob_intact:
                bucket_status[bucket_index] = 0
                num_perm += 1

        print num_perm
        timecourse[i] = num_perm
    rep_matrix[rep_index, :] = timecourse / num_buckets

rep_avgs = np.mean(rep_matrix, axis=0)
rep_stds = np.std(rep_matrix, axis=0)

time = np.arange(max_time)

plt.figure()
plt.errorbar(time, rep_avgs, yerr=rep_stds)

print rep_avgs[1]

"""
k = fitting.Parameter(0.1)
fmax = fitting.Parameter(0.5)
def exp_func(t):
    return fmax()*(1 - np.exp(-k()*t))
fitting.fit(exp_func, [k, fmax], timecourse, time)
plt.plot(time, exp_func(time))

m = fitting.Parameter(0.1)
def linear(t):
    return m()*t
fitting.fit(linear, [m], rep_avgs, time)
plt.plot(time, linear(time))

print m()
"""

