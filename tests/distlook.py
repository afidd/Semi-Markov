import logging
import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

logger=logging.getLogger(__file__)



def read_samples(fname):
    f=open(fname)
    x=list()
    y=list()
    for l in f.readlines():
        a,b=l.split()
        x.append(float(a))
        y.append(float(b))
    return np.array(x), np.array(y)


def test_dist(cdf_estimator):
    x, y=cdf_estimator
    #y2=scipy.stats.gamma.cdf(x, 1, scale=0.5)
    line, = plt.plot(x, y)
    #plt.plot(x,y2)
    plt.show()






if __name__=='__main__':
    samples=read_samples(sys.argv[1])
    test_dist(samples)
