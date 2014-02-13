import numpy as np


def linterp(m0, m1, p0, p1, p):
    """Interpolates linearly recarrays m0 and m1

    m0, m1 correspond to variables p0, p1
    Result is a recarray interpolated to the variable p (p0 < p < p1)
    """

    s = 1.*(p-p0)/(p1-p0)
    m = np.recarray(len(m0), dtype=m0.dtype)
    for col in m0.dtype.names:
        m[col] = (1.-s)*m0[col] + s*m1[col]
    return m

def linfit(x, y):
    """linear fit that returns only zero value, slope, and slope error
    """
    n = len(x)
    dp = n*sum(x**2) - (sum(x)**2)
    a = (sum(x**2)*sum(y) - sum(x)*sum(x*y))/dp
    b = (n*sum(x*y) - sum(x)*sum(y))/dp
    sigma = np.sqrt(sum((y-a-b*x)**2)/(n-2))
    err_a = np.sqrt((sigma**2/dp)*sum(x**2))
    err_b = np.sqrt((sigma**2/dp)*n)
    #return b, a, sigma, err_b, err_a
    return a, b, err_b
