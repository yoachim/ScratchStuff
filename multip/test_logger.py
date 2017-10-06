import numpy as np


class test_logger:
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


class supersimple:
    def __init__(self):
        self.dummy=5
    def __call__(self):
        return self.dummy
        