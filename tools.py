import numpy as np
import logging

logger = logging.getLogger(__name__)


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
    """Linear fit that returns only zero value, slope, and slope error
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

def read_csv(csv_file):
    """Reads CSV file with header and sends data to dictionary
    
    This routine was written specifically for q2 and is not meant 
    to be a generic CSV file reader. The CSV files that q2 reads must
    have an 'id' column. The first row must be a header column and it
    must not have empty cells. No commented rows are allowed. Empty
    data rows will be read as None. The output is a dictionary of numpy
    arrays.
    """
    f = open(csv_file)
    x = f.readlines()
    f.close()
    keys = x[0].split(",")
    if "" in keys:
        logger.error("First row of CSV file (keys) has an empty column.")
        return None
    if "id" not in keys:
        logger.error("CSV file must have an 'id' column.")
        return None
    dictionary = {}
    for key in keys:
        dictionary[key] = []
    for xi in x[1:]:
        for key, xij in zip(keys,xi.split(",")):
            if xij == "":
                xij = None
            if "teff" in key and xij != None:
                xij = int(xij)
            if ("logg" in key or "feh" in key or "vt" in key)\
               and xij != None:
                xij = float(xij)
            dictionary[key].append(xij)
    for key in keys:
        dictionary[key] = np.array(dictionary[key])
    ambiguous_ids = [idx for idx in dictionary["id"]\
                    if len(dictionary["id"][dictionary["id"]==idx]) > 1]
    if len(ambiguous_ids) > 1:
        logger.error("There are duplicates in 'id' column.")
        return None
    return dictionary
