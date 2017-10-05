import numpy as np
from multiprocessing import Pool


class test_logger(object):
    def __init__(self, i):
        self.i = i
        self.log = []

    def update_log(self, j):
        if np.max(j) > self.i:
            self.log.append(j)

    def __call__(self):
        if np.size(self.log) > 0:
            return np.mean(self.log)
        else:
            return 0


def update_log_helper(inobj, j):
    inobj.update_log(j)
    return inobj


def call_obj(inobj):
    return inobj()

if __name__ == '__main__':
    p = Pool(4)
    objs = [test_logger(i) for i in range(4)]

    numbers = [2, 2, 4, 5, 2, 3, 5, 6, 1, 1, 1]
    for number in numbers:
        l = [np.array([number,number+1])]*len(objs)
        objs = p.starmap(update_log_helper, zip(objs, l))

    results = p.map(call_obj, objs)
    print(results)
