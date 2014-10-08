#/opt/local/bin/python3.3

#
# Formulae taken fro Ball, F. and Nasell, I. 'The shape of the size 
# distribution of an epidemic in a finite population.' 
# Math. Biosc. 123:167-181 (1994).
#

import sys, argparse
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import matplotlib.pylab as pylab
import collections

def indices(N, a):
    indx = collections.OrderedDict()
    n = 0
    for i in range(N+1):
        for j in range(1, N+a-i+1):
            assert(0 <= i and i <= N and 1 <= i+j and i+j <= N+a)
            indx[(i, j)] = n
            n += 1
    return indx

def Bailey(N, a, R0, indx):
    row = []
    col = []
    data = []

    m = indx[(N, a)]
    row.append(m)
    col.append(m)
    data.append(a*N*(1.0+1.0/R0))

    for (i, j), m in indx.items():
        if i+1 <= N and 2 <= j:
            row.append(m)
            col.append(indx[(i+1,j-1)])
            data.append((i+1.0)*(j-1.0))
        
        if i != N or j != a:
            row.append(m)
            col.append(m)
            data.append(-j*(i+N/R0))

        if i+j+1 <= N+a:
            row.append(m)
            col.append(indx[(i,j+1)])
            data.append(N*(j+1.0)/R0)
           
    return sparse.csr_matrix(sparse.coo_matrix((data, (row, col))))

def rhs(N, a, indx):
    v = []
    for (i, j), n in indx.items():
        if i == N and j == a:
            t = 1.0
        else:
            t = 0.0
        v.insert(n, t)

    return v

class CustomParser(argparse.ArgumentParser):
    
    def __init__(self):
        super(CustomParser, self).__init__(description='Bailey\'s final size distribution for SIR model',
                                           epilog='Formulae taken from Ball, F. and Nasell, I. "The shape of the size distribution of an epidemic in a finite population." Math. Biosc. 123:167-181 (1994).')

        self.add_argument('-N', '--susceptible-population', type=int, default=49,
                           dest='N', metavar='<count>', 
                           help='Initial number of susceptible hosts')
        self.add_argument('-a', '--infected-population', type=int, default=1,
                          dest='a', metavar='<count>', 
                          help='Initial number of infected hosts')
        self.add_argument('-r', '--basic-reproduction-number', type=float, default=3,
                          dest='R0', metavar='<ratio>',
                          help='Basic reproduction number')
        self.add_argument('-s', '--sparse-matrix-plot', dest='splot', action='store_true', default=False)
        self.add_argument('-p', '--probability-plot', dest='pplot', action='store_true', default=False)
        self.add_argument('-i', '--list-indices', dest='ilist', action='store_true', default=False)
        self.add_argument('-m', '--mean', dest='mean', action='store_true', default=False)
        
if __name__ == '__main__':
    parser = CustomParser()
    ns = parser.parse_args(sys.argv[1:])

    indx = indices(ns.N, ns.a)

    if ns.ilist:
        for (i, j), n in indx.items():
            print(i, j, n)

    B = Bailey(ns.N, ns.a, ns.R0, indx)
    d = rhs(ns.N, ns.a, indx)

    f = linalg.spsolve(B, np.array(d))

    x = []
    y = []
    for (i, j), n in indx.items():
        if j==1:
            x.append(ns.N-i)
            y.append(ns.N*f[n]/ns.R0)

    assert(abs(1.0-sum(y)) < 1.0e-12)

    if ns.mean:
        print('sum %f' % (np.sum(np.array(y)),))
        print('mean %f' % (np.sum((ns.a+np.array(x))*np.array(y)),))

    if ns.splot:
        pylab.figure()
        pylab.spy(B)
        pylab.title('Sparsity pattern')

    if ns.pplot:
        pylab.figure()
        pylab.plot(x, y)
        pylab.title('Probability mass function of outbreak sizes\nN=%d, a=%d, R0=%f (dimension %d)' % (ns.N, ns.a, ns.R0, len(indx)))

    print('# N=%d, a=%d, R0=%f, total dimension of matrix is %d' % (ns.N, ns.a, ns.R0, len(indx)))
    for i in range(len(x)):
        print('%d, %f' % (int(x[i]), y[i]))

    if ns.splot or ns.pplot:
        pylab.show()

    sys.exit(0)

