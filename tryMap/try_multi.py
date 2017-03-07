import numpy as np
from multiprocessing import Pool


def run_con(time_series_ker):
    time_series_ker[0].convolve(time_series_ker[1])
    return time_series_ker[0]


class time_series(object):
    def __init__(self, values):
        self.t = values

    def convolve(self, kernel=10):
        kernel = np.ones(kernel+2)
        kernel[0] = 0
        kernel[-1] = 0
        self.smooth = np.convolve(self.t, kernel)

if __name__ == "__main__":
    lengths = range(100, 105)
    tss = [time_series(np.random.rand(length)) for length in lengths]
    pool = Pool(3)
    result = pool.map(run_con, zip(tss, lengths))
    