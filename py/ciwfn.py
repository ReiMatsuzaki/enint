import pandas as pd
import numpy as np
tr = np.transpose
dot = np.dot

from enint.math import ijv2mat

class Aij(object):
    """
    manage one particle coupling constant
    a_{IJ}^{ij} = <Phi_I|E_{ij}|Phi_J>
    """
    def __init__(self, ilist, jlist, iilist, jjlist, vlist, nfrozen):
        self.ilist = list(ilist)
        self.jlist = list(jlist)
        self.iilist = list(iilist)
        self.jjlist = list(jjlist)
        self.vlist = list(vlist)
        self.ncsf = max(iilist)

        for i in range(1,nfrozen+1):
            for ii in range(1,self.ncsf+1):
                self.ilist.append(i)
                self.jlist.append(i)
                self.iilist.append(ii)
                self.jjlist.append(ii)
                self.vlist.append(2.0)
                
        self.ilist  = np.array(self.ilist, dtype=int)
        self.jlist  = np.array(self.jlist, dtype=int)
        self.iilist = np.array(self.iilist, dtype=int)
        self.jjlist = np.array(self.jjlist, dtype=int)
        self.vlist  = np.array(self.vlist)


    def mo2csf(self, mmo):
        mcsf=np.zeros((self.ncsf,self.ncsf))
        for (i,j,ii,jj,v) in zip(self.ilist,
                                 self.jlist,
                                 self.iilist,
                                 self.jjlist,
                                 self.vlist):
            if(ii==jj and i==j):
                mcsf[ii-1,jj-1] += v*mmo[i-1,j-1]
            elif(ii!=jj and i!=j):
                mcsf[ii-1,jj-1] += v*mmo[i-1,j-1]
                mcsf[jj-1,ii-1] += v*mmo[j-1,i-1]
            else:
                raise RuntimeError("illegal combination of i,j,I,J")
        return mcsf

    def mo2ci(self, mmo, cci):
        nstate = len(cci[0,:])
        mci = np.zeros((nstate, nstate))
        for (i,j,ii,jj,v) in zip(self.ilist,
                                 self.jlist,
                                 self.iilist,
                                 self.jjlist,
                                 self.vlist):
            if(ii==jj and i==j):
                for n in range(nstate):
                    for m in range(nstate):
                        mci[n,m] += cci[ii-1,n]*cci[jj-1,m]*v*mmo[i-1,j-1]
            elif(ii!=jj):
                for n in range(nstate):
                    for m in range(nstate):
                        mci[n,m] += cci[ii-1,n]*cci[jj-1,m]*v*mmo[i-1,j-1]
                        mci[n,m] += cci[jj-1,n]*cci[ii-1,m]*v*mmo[j-1,i-1]
            else:
                raise RuntimeError("illegal combination of i,j,I,J")
        return mci
        

def aij_load(fn, nfrozen):
    df = pd.read_csv(fn)
    return Aij(df["i"], df["j"], df["I"], df["J"], df["val"], nfrozen)

def ao2mo(mao, cmo):
    """
    Inputs
    ------
    mao : matrix
    .     AO represented matrix
    cmo : matrix
    .     MO coefficient matrix

    Returns
    -------
    mmo : matrix
    .     MO represented matrix
    """
    return dot(tr(cmo), dot(mao, cmo))




