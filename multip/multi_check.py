import numpy as np
from multiprocessing import Pool


class test_logger(object):
    def __init__(self, maxval):
        self.maxval = maxval
        self.log = []

    def update_log(self, j):
        if np.max(j) > self.maxval:
            self.log.append(j)

    def __call__(self):
        if np.size(self.log) > 0:
            return np.std(self.log)
        else:
            return 0


def update_log_helper(inobj, j):
    inobj.update_log(j)
    return inobj


def call_obj(inobj):
    return inobj()

if __name__ == '__main__':
    p = Pool(4)
    loggers = [test_logger(i) for i in range(4)]

    numbers = [2, 2, 4, 5, 2, 3, 5, 6, 1, 1, 1]
    for number in numbers:
        l = [np.array([number, number+1])]*len(loggers)
        loggers = p.starmap(update_log_helper, zip(loggers, l))

    results = p.map(call_obj, loggers)
    print(results)
