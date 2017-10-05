import numpy as np
from multiprocessing import Pool


## Try out multiprocessing with some objects


class test_object(object):
    def __init__(self, i):
        self.i = i

    def __call__(self, j):
        return self.i * j

def run_object(inobj):
    return inobj(5.)


if __name__ == '__main__':
    p = Pool(4)
    objs = [test_object(i) for i in range(4)]

    results = p.map(run_object, objs)
    print(results)
