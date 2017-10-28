import types

import numpy as np
import pandas as pd

dot = np.dot
def ijv2mat(df):
    ilist = df["i"]
    jlist = df["j"]
    n = max(ilist)
    m = max(jlist)
    
    if("val" in df):
        vlist = df["val"]
        dtype = float
    elif("re" in df and "im" in df):
        vlist = df["re"] + 1.0j * df["im"]
        dtype = complex

    mat = np.zeros((n,m), dtype=dtype)
    for (i,j,v) in zip(ilist,jlist,vlist):
        mat[i-1,j-1]=v
    return mat

def mat2ijv(x, eps=None):
    (n,m) = np.shape(x)
    ilist = []
    jlist = []
    vlist = []
    idx = 0
    for i in range(n):
        for j in range(m):
            if(eps==None):
                non0 = True
            else:
                non0 = (abs(x[i,j])>eps)
            if(non0):
                ilist.append(i+1)
                jlist.append(j+1)
                vlist.append(x[i,j])
                
    i = np.array(ilist)
    j = np.array(jlist)
    v = np.array(vlist)    
    if(isinstance(x[0,0], complex)):
        df = pd.DataFrame({"i":i, "j":j, "re":v.real, "im":v.imag},
                          columns=["i","j","re","im"])
    else:
        df = pd.DataFrame({"i":i, "j":j, "val":v},
                          columns=["i","j","val"])        
    return df

def iv2vec(df):
    ilist = df["i"]
    n = max(ilist)
    if("val" in df):        
        vlist = df["val"]
        vec = np.zeros(n,dtype=float)
    elif("re" in df and "im" in df):        
        vlist = df["re"] + 1.0j * df["im"]
        vec = np.zeros(n,dtype=complex)
    for (i,v) in zip(ilist, vlist):
        vec[i-1] = v
    return vec

def vec2iv(vec):
    n = len(vec)
    ilist = []
    vlist = []
    for i in range(n):
        ilist.append(i+1)
        vlist.append(vec[i])
    if(isinstance(vec[0], complex)):
        val = np.array(vlist)
        re = val.real
        im = val.imag
        i = np.array(ilist)
        df = pd.DataFrame({"i":i, "re":re, "im":im},
                          columns=["i","re","im"])
    else:
        val = np.array(vlist)
        df = pd.DataFrame({ "i": np.array(ilist), "val": val})
    return df

if __name__ == '__main__':
    x = np.array([[2.0, 1.1],
                  [2.1, 1.2],
                  [3.2, 1.3]])
    df = mat2ijv(x)
    y = ijv2mat(df)
    print sum(abs(y-x))


def hdot(x, y):
    return dot(x.conj(), y)

def norm(x):
    return np.sqrt(hdot(x,x).real)
    
