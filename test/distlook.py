import logging
import sys
import numpy as np
import scipy.stats
import matplotlib
print(matplotlib.get_backend())
matplotlib.use('pdf')
import matplotlib.pyplot as plt

logger=logging.getLogger(__file__)



def read_samples(fname):
    f=open(fname)
    x=list()
    y=list()
    for l in f.readlines():
        a,b=l.split()
        try:
            x.append(float(a))
            y.append(float(b))
        except:
            pass
    return np.array(x), np.array(y)


def test_dist(cdf_estimator):
    x, y=cdf_estimator
    #y2=scipy.stats.gamma.cdf(x, 1, scale=0.5)
    line, = plt.plot(x, y)
    #plt.plot(x,y2)
    #plt.draw()

    plt.show()


def compare_dist(a, b):
    logger.debug('{0} samples'.format(len(a[0])))
    x, y=a
    #y2=scipy.stats.gamma.cdf(x, 1, scale=0.5)
    line, = plt.plot(x, y)
    x2, y2=b
    plt.plot(x2, y2, '-')
    plt.draw()
    plt.savefig("fig.pdf")
    plt.show()



if __name__=='__main__':
    logging.basicConfig(level=logging.DEBUG)
    samples=read_samples(sys.argv[1])
    if len(sys.argv)==2:
        test_dist(samples)
    elif len(sys.argv)==3:
        logger.info("Comparing two distributions")
        #python3 tests/distlook.py piecewise.txt stdpiece.txt
        compare=read_samples(sys.argv[2])
        compare_dist(samples, compare)
    else:
        print("Looking for one or two arguments.")
