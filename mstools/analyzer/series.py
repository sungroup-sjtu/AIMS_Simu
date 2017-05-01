import random

from pandas import Series
import numpy as np


def is_converged(series: Series, tolerance=0.9) -> (bool, float):
    n_block = 5
    start = series.index[0]
    end = series.index[-1]
    span = end - start
    block_size = span / n_block

    blocks = []
    dists = []
    for i in range(n_block):
        block = series[start + i * block_size: end + i * block_size]
        blocks.append(block)
        mu = np.mean(block)
        sigma = np.std(block)
        dists.append([mu, sigma])

    minimal = dists[0][0] - dists[0][1] * 3.8
    maximal = dists[0][0] + dists[0][1] * 3.8
    height = max([0.3989 / dist[1] for dist in dists])
    for dist in dists[1:]:
        minimal = min(minimal, dist[0] - dist[1] * 3.8)
        maximal = max(maximal, dist[0] + dist[1] * 3.8)

    points = int(1E5)
    random_x_list = np.ndarray(shape=points)
    random_y_list = np.ndarray(shape=points)
    gauss_y_list = np.ndarray(shape=(points, n_block), dtype=float)
    for n in range(points):
        x = random.random() * (maximal - minimal) + minimal
        y = random.random() * height
        random_x_list[n] = x
        random_y_list[n] = y
        for i in range(n_block):
            gauss_y_list[n][i] = gauss(dists[i][0], dists[i][1], x)

    # Plot the points
    # import pylab
    # from scipy.stats import norm
    # pylab.xlim([minimal, maximal])
    # pylab.ylim([0, height])
    # pylab.scatter(random_x_list, random_y_list, s=1)
    # for i in range(n_block):
    #     dist = dists[i]
    #     x_series = np.linspace(minimal, maximal, 100)
    #     pylab.plot(x_series, norm.pdf(x_series, dist[0], dist[1]), linewidth=1 + i * 0.5)
    # pylab.show()

    for i in range(n_block - 2):
        overlap = 0
        for n in range(points):
            if random_y_list[n] < min(gauss_y_list[n][i:]):
                overlap += 1
        overlap /= points
        overlap *= ((maximal - minimal) * height)
        print(i, overlap)
        if overlap > tolerance:
            converged_from = start + block_size * i
            return True, converged_from

    else:
        return False, 0


def gauss(mu, sigma, x):
    h = 0.3989 / sigma
    exponent = -(x - mu) ** 2 / 2 / sigma ** 2
    return h * np.exp(exponent)

def average(series)->(float, float):
    '''
    get average and standard error
    '''
    pass

