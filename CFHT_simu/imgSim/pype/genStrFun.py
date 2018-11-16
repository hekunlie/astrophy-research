# General Functions for String Operations.
# Liu Dezi

from astropy.io import fits
import os, sys

def read_param(filename):
    """
    Read a parameter configuration file, and save the parameters into
    a python directory structure.

    Parameter:
    filename: the input configuration file name.

    NOTE: The format of the configuration should be:
    -----------------------------------------------
    ########### Configuration file ################
    # blank or some comments
    param1 value # comment
    param2 value
    ...

    # blank or some comments
    paramn value
    ...
    ------------------------------------------------
    For each parameter, only one value (string, float, integer ...) can
    be assigned. You have to avoid double '#' (i.e. '##'), but it will
    work for many more '#' or a single '#'.
    """

    pfile = open(filename).readlines()
    nn = len(pfile)
    param = {} # define dictionary structure
    for i in range(nn):
        rowlist = pfile[i].split()
        if len(rowlist)<=1: continue # blank row
        if not "#" in rowlist:
            if len(rowlist)==2:
                key, value = rowlist[0:2]
                param.update({key:value})
            else:
                print "!! Something is wrong with parameter '%s'."%rowlist[0]
                return
        elif rowlist.index("#")==2:
            key, value = rowlist[0:2]
            param.update({key:value})
        elif rowlist.index("#")==0:
            continue # annotation
        else:
            print "!! Something is wrong with parameter '%s'."%rowlist[0]
            return
    return param

def write_FitsTable(table_contents, table_array, outcat="table.ldac"):
    """
    Write FITS table
    """
    conts = open(table_contents).read().splitlines()
    rid = [i for i in range(len(conts)) if not "#" in conts[i]]
    ncont = len(rid)

    x, y = table_array.shape
    if y!=ncont: raise ValueError("!!! Dimensition error")

    # first construct an empty fits table
    prihdr = fits.Header()
    prihdr["AUTHOR"]="Dezi LIU, Peking University"
    prihdu = fits.PrimaryHDU(header=prihdr)

    # Table Extension
    cols = []
    for i in range(ncont):
        cont = conts[rid[i]]
        cont = cont.split()
        name, fmt, disp = cont[0], cont[1], " ".join(cont[2:])
        col = fits.Column(name=name, format=fmt, disp=disp, array=table_array[:,i])
        cols += [col]

    cols = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.header["EXTNAME"] = "OBJECTS"

    thdulist = fits.HDUList([prihdu,tbhdu])
    if os.path.exists(outcat): os.popen("rm %s"%outcat)
    thdulist.writeto(outcat)

    return

