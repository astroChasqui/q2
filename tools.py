import numpy as np
import logging
import os

logger = logging.getLogger(__name__)


def linterp(m0, m1, p0, p1, p):
    """Interpolates linearly dictionaries of numpy arrays m0 and m1

    m0, m1 correspond to variables p0, p1
    Result is a dictionary interpolated to the variable p (p0 < p < p1)
    """
    s = 1.*(p-p0)/(p1-p0)
    m = {}
    for col in m0.keys():
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

def read_csv(csv_file, file_type=None):
    """Reads CSV file with header and sends data to dictionary
    
    file_type = None, 'stars', 'lines'

    This routine was written specifically for q2 and is not meant 
    to be a generic CSV file reader. The first row must be a header column
    and it must not have empty cells. No commented rows are allowed. Empty
    data rows will be read as None. The output is a dictionary of numpy
    arrays.
    
    If file_type is None, it does not matter what the header contains,
    however, all data must be numbers.
    If file_type is 'stars', the CSV file must have an 'id' column.
    If file_type is 'lines', the CSV file must have at least these colums:
    'wavelength', 'species', 'ep', and 'gf'.
    """

    with open(csv_file, 'rU') as f:
        x = f.readlines()

    keys = [key.strip("\n") for key in x[0].split(",")]
    if (len(keys) != len(set(keys))):
        logger.error("First row of CSV file (keys in "+csv_file+\
                     ") has duplicates.")
        return None
    if "" in keys:
        logger.error("First row of CSV file (keys in "+csv_file+\
                     ") has empty columns.")
        return None

    if file_type == "stars":
        if "id" not in keys:
            logger.error("Stars CSV file must have an 'id' column.")
            return None
        important_pars = ["teff", "logg", "feh"]
        for par in important_pars:
            if par not in keys:
                logger.warning("Stars CSV file does not have a '"+par+\
                               "' column.")

    if ("wavelength" not in keys or "species" not in keys\
        or "ep" not in keys or "gf" not in keys) and file_type=="lines":
        logger.error("Lines CSV file must have all of these columns: "+\
                     "'wavelength', 'species', 'ep', and 'gf'.")
        return None

    dictionary = {}
    for key in keys:
        dictionary[key] = []
    for xi in x[1:]:
        for key, xij in zip(keys,xi.split(",")):
            xij = xij.strip("\n")
            if xij == "":
                xij = None
            if "teff" in key and xij != None:
                xij = int(xij)
            if ("logg" in key or "feh" in key or "vt" in key)\
               and xij != None:
                xij = float(xij)
            if (key == "v" or "err" in key or "plx" in key)\
               and xij != None:
                xij = float(xij)
            if file_type == "lines" and xij != None and key != "comments":
                xij = float(xij)
            if file_type == None:
                xij = float(xij)
            dictionary[key].append(xij)
    for key in keys:
        dictionary[key] = np.array(dictionary[key])
    if file_type == "stars":
        if [None] in dictionary["id"]:
            logger.error("The 'id' column cannot have empty rows.")
            return None
        ambiguous_ids = [idx for idx in dictionary["id"]\
                        if len(dictionary["id"][dictionary["id"]==idx]) > 1]
        if len(ambiguous_ids) > 1:
            logger.error("There are duplicates in 'id' column.")
            return None
    return dictionary
