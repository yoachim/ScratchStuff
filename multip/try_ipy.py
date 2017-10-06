import numpy as np
import ipyparallel as ipp

# OK, let's try this with the ipython parallel stuff

# ipcluster start -n 4



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
    rc = ipp.Client()

    dview = rc[:]
    loggers = [test_logger(i) for i in range(4)]
    # put one logger on each engine
    for i, logger in enumerate(loggers):
        rc[i].push({'logger':logger})

    numbers = [2, 2, 4, 5, 2, 3, 5, 6, 1, 1, 1]
    # Let's update the objects
    
    import pdb ; pdb.set_trace()
    for number in numbers:
        dview.apply_sync(lambda x: logger.update_log(x), number)


    dview.execute('result = logger()')
    results = dview['result']
    print(results)
    

    #for number in numbers:
    #    l = [np.array([number, number+1])]*len(loggers)
    #    loggers = p.starmap(update_log_helper, zip(loggers, l))

    #results = p.map(call_obj, loggers)
    #print(results)
