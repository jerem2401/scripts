#!/usr/bin/env python
#script given by Giovanni Bussi#

import re
import gzip

class FormatError(Exception):
    """Custom error reported by read_as_pandas.
    """
    pass

def read_as_pandas(file,chunksize=None,usecols=None,skiprows=None,nrows=None):
    """Import a plumed data file as a pandas dataset.

       file : either path to a file or open filed object

       chunksize : int, optional
           Return an iterable object.
       usecols : list-like or callable, optional
           Directly passed to pandas.
       skiprows : list-like, int or callable, optional
           Directly passed to pandas.
       nrows : int, optional
           Directly passed to pandas.

       Returns
       -------
       DataFrame (when chunksize is not provided) or iterable TextFileReader (when chunksize is provided).

       Comments
       --------

       Gzipped files are supported and automatically detected when a file name ends with '.gz'.

       `pandas` module is imported the first time this function is used. Since importing `pandas` is quite slow,
       the first call to this function will be significantly slower than the following ones.
       Following calls should be fast. The overall speed is comparable or better to loading with `numpy.loadtxt`.

       Examples
       --------

       colvar=plumed.read_as_pandas("COLVAR")
       print(colvar) # show the datasheet

       colvar=plumed.read_as_pandas("COLVAR",usecols=[0,4])
       print(colvar) # show the datasheet

       colvar=plumed.read_as_pandas("COLVAR",usecols=["time","distance"])
       print(colvar) # show the datasheet

       for chunk in plumed.read_as_pandas("COLVAR",chunksize=10):
           print(chunk) # show the datasheet. actually here you should process the chunk

       Limitations
       -----------

       1. Constants are presently ignored. As a consequence it is not possible to retrieve
       information such as ranges of variables.

       2. Text variables are not converted using standard PLUMED conversions (e.g. `pi` converted to
       the value of pi and arithmetic operations resolved).

       3. Only the initial header is read, which implies that files resulting from concatenating
       datasets with a different number of columns or different columnd names will not
       be read correctly.

       Issue 1 could be solved parsing the `#! SET` lines in python.
       Issue 2 could be solved calling `PLMD::Tools::convert`, which
       could be easily done through the `plumed_cmd` interface (TODO).
 
       Alternatively, all issues might be solved using `PLMD::IFile` for reading,
       which could be useful but possibly more complicated to implement.
    """
# importing pandas is pretty slow, so we only do it when needed
    import pandas as pd
# allow passing a string
    if isinstance(file,str):
        file=open(file)
# take care of gzipped files
    if re.match(".*\.gz",file.name):
        file = gzip.open(file.name,'rt')
# read first line
    line = file.readline()
    columns = line.split()
# check header
    if len(columns)<2:
        raise FormatError("Error reading PLUMED file "+file.name + ". Not enough columns")
    if columns[0] != "#!" or columns[1] != "FIELDS":
        raise FormatError("Error reading PLUMED file" +file.name + ". Columns: "+columns[0]+" "+columns[1])
# read column names
    columns = columns[2:]
# read the rest of the file
# notice that if chunksize was provided the result will be an iterable object
    return pd.read_csv(file, delim_whitespace=True, comment="#", header=None,names=columns,
                       usecols=usecols,skiprows=skiprows,nrows=nrows,chunksize=chunksize)
