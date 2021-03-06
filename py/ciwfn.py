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
    def __init__(self, ilist, jlist, iilist, jjlist, vlist):
        self.ilist = list(ilist)
        self.jlist = list(jlist)
        self.iilist = list(iilist)
        self.jjlist = list(jjlist)
        self.vlist = list(vlist)
        
        self.ncsf = max(iilist)
        self.nmo = max(ilist)

        self.nfrozen = min(ilist)-1

        """
        for i in range(1,self.nfrozen+1):
            for ii in range(1,self.ncsf+1):
                self.ilist.append(i)
                self.jlist.append(i)
                self.iilist.append(ii)
                self.jjlist.append(ii)
                self.vlist.append(2.0)
        """
                
        self.ilist  = np.array(self.ilist, dtype=int)
        self.jlist  = np.array(self.jlist, dtype=int)
        self.iilist = np.array(self.iilist, dtype=int)
        self.jjlist = np.array(self.jjlist, dtype=int)
        self.vlist  = np.array(self.vlist)

    def mo2csf(self, mmo):
        mcsf=np.zeros((self.ncsf,self.ncsf))

        for i in range(1,self.nfrozen+1):
            for ii in range(1,self.ncsf+1):
                mcsf[ii-1,ii-1] += 2*mmo[i-1,i-1]
        
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

        for i in range(1,self.nfrozen+1):
            for ii in range(1,self.ncsf+1):
                for n in range(nstate):
                    for m in range(nstate):
                        mci[n,m] += cci[ii-1,n]*cci[ii-1,m]*2*mmo[i-1,i-1]
        
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

    def dm1(self, in_bcci, in_kcci):
        """
        compute one particle density matrix
        Psi = \sum cI Phi_I  
        where cI is CI coefficient

        Inputs
        ------
        cci : Matrix
        .    cci[:] means n th CI vector

        Returns
        -------
        dm1 : ndarray(nmo,nmo,nci)
        .    dm1n[:,:] means n th density matrix
        """        

        if(isinstance(in_kcci[0], float)):
            kcci = in_kcci
            bcci = in_bcci
            dm1 = np.zeros((self.nmo, self.nmo))
        elif(isinstance(in_kcci[0], complex)):
            kcci = in_kcci
            bcci = in_bcci.conj()
            dm1 = np.zeros((self.nmo, self.nmo), dtype=complex)
        else:
            raise RuntimeError("invalid cci")

        for i in range(1,self.nfrozen+1):
                for ii in range(1,self.ncsf+1):
                    dm1[i-1,i-1] += bcci[ii-1]*kcci[ii-1]*2
            
        for (i,j,ii,jj,aijIJ) in zip(self.ilist,
                                 self.jlist,
                                 self.iilist,
                                 self.jjlist,
                                 self.vlist):
            if(ii==jj and i==j):
                dm1[i-1,j-1] += bcci[ii-1]*kcci[jj-1]*aijIJ
            elif(ii!=jj and i!=j):
                dm1[i-1,j-1] += bcci[ii-1]*kcci[jj-1]*aijIJ
                dm1[j-1,i-1] += bcci[ii-1]*kcci[jj-1]*aijIJ
            else:
                raise RuntimeError("illegal combination of i,j,I,J")
        return dm1

def ciwfn_op1(mmo, dm1):
    """
    compute expectation value of one particle operator
    <Psi_n | O | Psi_n> = sum_{ij} o_ij D^n_ij
    """
    nmo = len(dm1[:,0])
    return np.sum(mmo[:nmo,:nmo]*dm1[:,:])
        
def aij_load(fn):
    df = pd.read_csv(fn)
    return Aij(df["i"], df["j"], df["I"], df["J"], df["val"])

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

def rho_mo(psi_ao, cmo):
    psi_mo = np.dot(tr(cmo), psi_ao) 
    rho_mo = psi_mo[np.newaxis, :, :] * psi_mo[:, np.newaxis, :]
    return rho_mo

    
