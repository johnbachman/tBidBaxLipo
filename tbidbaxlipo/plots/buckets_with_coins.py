import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from tbidbaxlipo.util import fitting
import sys

plt.ion()

def exp_fits():
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

    # Do a series of replicates to get meaningful averages
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

def scaling_of_average():
    """Goal is to look at scaling of average coins per bucket when the
    distribution of coins among buckets is Poisson vs. when it is more likely
    that coins go to buckets that already have coins.

    New algorithm:

    Each coin chooses a bucket from 1 to n. After choosing a bucket, you
    flip a coin to see if the coin goes to that bucket. The bias of the
    coin increases with the number of coins already in that bucket. If the
    coin is heads, you assign the coin to that bucket.

    Old algorithm:

    For each coin, a probability that it stays undistributed or it goes to a
    bucket. If it goes to a bucket, it is randomly assigned to a bucket.  If no
    coins in any buckets, you choose a number from 1 to n (buckets).  Imagine a
    list of 1...100. Now, imagine that bucket number 3 has a coin.  Suppose it
    is now twice as likely to go that bucket: you could put 1, 2, 3, 3, 4,
    ...100. In the bucket list so that 3 is twice as likely to be chosen. Then,
    for every additional coin, you update the bucket list accordingly. This can
    be done by simply appending to the end.

    Algorithm stops once all coins have been distributed. At the end, you count
    the number of entries in the bucket list, and subtract 1 from all the counts.
    This gives you the number of coins in that bucket.
    """
    coin_arr = np.arange(0, 600, 4)
    num_avgs = coin_arr.shape[0]
    num_rounds = 10
    averages = np.zeros((num_avgs, num_rounds))
    num_buckets = 100.

    for i, num_coins in enumerate(coin_arr):
        for j in range(num_rounds):
            coin_status = np.zeros(num_coins, dtype='int')
            p_bind = 0.2
            while np.count_nonzero(coin_status) < num_coins:
                for coin_index in np.where(coin_status == 0)[0]:
                    bucket_index = np.random.randint(1, num_buckets+1)
                    p = np.random.rand()
                    if p < p_bind:
                        coin_status[coin_index] = bucket_index
            bucket_counts = np.bincount(coin_status)
            avg_nonzero = np.mean(bucket_counts[np.nonzero(bucket_counts)])
            averages[i, j] = avg_nonzero
    plt.figure()
    plt.errorbar(coin_arr,
                 np.mean(averages, axis=1),
                 yerr=np.std(averages, axis=1) / np.sqrt(num_rounds))
    plt.plot(coin_arr, coin_arr/num_buckets)
    plt.xlabel('Number of total coins')
    plt.ylabel('Average coins in non-empty buckets')
    plt.title('Poisson distribution')
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    scaling_of_average()

