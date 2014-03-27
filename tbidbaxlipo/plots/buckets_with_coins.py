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
                # If our random number is greater than prob intact, we
                # permeabilize
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

def scaling_of_average(model='poisson'):
    """Goal is to look at scaling of average coins per bucket when the
    distribution of coins among buckets is Poisson vs. when it is more likely
    that coins go to buckets that already have coins.

    New algorithm:

    Each coin chooses a bucket from 1 to n. After choosing a bucket, you
    flip a coin to see if the coin goes to that bucket. The bias of the
    coin increases with the number of coins already in that bucket. If the
    coin is heads, you assign the coin to that bucket.

    Algorithm stops once all coins have been distributed. At the end, you count
    the number of coins in the buckets.
    """
    # The different amounts of total coins that we will try
    total_coin_arr = np.arange(0, 2000, 10)
    # The number of different experiments with different total coins
    num_coin_totals = total_coin_arr.shape[0]
    # Number of rounds (replicates) for each coin number to establish estimate
    num_rounds = 10
    # Array for storing the calculated averages
    averages = np.zeros((num_coin_totals, num_rounds))
    # The number of buckets into which we distribute coins
    num_buckets = 100

    # Try each coin total
    for i, num_coins in enumerate(total_coin_arr):
        # For each coin total do num_rounds replicates
        for j in range(num_rounds):
            coin_status = np.zeros(num_coins)
            coin_status[:] = np.nan
            # Baseline probability of binding
            p_bind_base = 0.1
            p_nobind_base = 1 - p_bind_base
            bucket_counts = np.zeros(num_buckets)
            while np.isnan(coin_status).sum() > 0:
                for coin_index in np.where(np.isnan(coin_status))[0]:
                    bucket_index = np.random.randint(num_buckets)
                    num_coins_in_bucket = bucket_counts[bucket_index]
                    # Choose our model, either autoactivation or
                    # poisson (default).
                    # Autoactivation model: we make it such that for every
                    # additional coins that there is in the bucket, the
                    # probability of NOT binding goes down by half (relative to
                    # baseline). So the first coin will obviously make a much
                    # bigger difference in absolute terms than the 10th coin.
                    if model == 'auto':
                        p_bind = 1 - (p_nobind_base *
                                      (0.5) ** num_coins_in_bucket)
                    else:
                        p_bind = p_bind_base
                    # Decide whether the coin "inserts" or not
                    p = np.random.rand()
                    if p < p_bind:
                        coin_status[coin_index] = bucket_index
                        bucket_counts[bucket_index] += 1
            # Calculate the average of the nonzero entries
            avg_nonzero = np.mean(bucket_counts[np.nonzero(bucket_counts)])
            averages[i, j] = avg_nonzero
    # Plot average
    plt.figure()
    plt.errorbar(total_coin_arr,
                 np.mean(averages, axis=1),
                 yerr=np.std(averages, axis=1) / np.sqrt(num_rounds))
    plt.plot(total_coin_arr, total_coin_arr/float(num_buckets))
    plt.xlabel('Number of total coins')
    plt.ylabel('Average coins in non-empty buckets')
    plt.title('Poisson distribution')
    import ipdb; ipdb.set_trace()

if __name__ == '__main__':
    scaling_of_average(model='auto')

