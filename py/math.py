import types

import numpy as np
import pandas as pd

dot = np.dot
def gauss(a,zeta,r0,p0,x):
    return a*np.exp(-zeta*(x-r0)**2 + 1.0j*p0*(x-r0))

def eigh_sort(h):
    (e,u) = np.linalg.eigh(h)
    eu_list = []
    for i in range(len(e)):
        eu_list.append((e[i], u[:,i]))
    eu_list.sort(key=lambda eu:eu[0])
    e = np.array([ei for (ei,ui) in eu_list])
    u = np.transpose(np.array([ui for (ei,ui) in eu_list]))
    return (e,u)

def ijkv2ten(df):
    
    if(isinstance(df, str)):
        df2 = pd.read_csv(df)
        return ijv2mat(df2)
    
    ilist = df["i"]
    jlist = df["j"]
    klist = df["j"]
    ni = max(ilist)
    nj = max(jlist)
    nk = max(klist)
    
    if("val" in df):
        vlist = df["val"]
        dtype = float
    elif("re" in df and "im" in df):
        vlist = df["re"] + 1.0j * df["im"]
        dtype = complex

    ten = np.zeros((ni,nj,nk), dtype=dtype)
    for (i,j,k,v) in zip(ilist,jlist,klist,vlist):
        ten[i-1,j-1,k-1]=v
    return ten

def ten2ijkv(x, eps=None):
    (ni,nj,nk) = np.shape(x)
    ilist = []
    jlist = []
    klist = []
    vlist = []
    idx = 0
    for i in range(ni):
        for j in range(nj):
            for k in range(nk):
                if(eps==None):
                    non0 = True
                else:
                    non0 = (abs(x[i,j,k])>eps)
                if(non0):
                    ilist.append(i+1)
                    jlist.append(j+1)
                    klist.append(k+1)
                    vlist.append(x[i,j,k])
                
    i = np.array(ilist)
    j = np.array(jlist)
    k = np.array(klist)    
    v = np.array(vlist)    
    if(isinstance(x[0,0,0], complex)):
        df = pd.DataFrame({"i":i, "j":j, "k":k, "re":v.real, "im":v.imag},
                          columns=["i","j","k","re","im"])
    else:
        df = pd.DataFrame({"i":i, "j":j, "k":k, "val":v},
                          columns=["i","j","k","val"])        
    return df

def ijv2mat(df):

    if(isinstance(df, str)):
        df2 = pd.read_csv(df)
        return ijv2mat(df2)
    
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

    if(isinstance(df, str)):
        df2 = pd.read_csv(df)
        return iv2vec(df2)
    
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
    print(sum(abs(y-x)))



def hdot(x, y):
    return dot(x.conj(), y)

def norm(x):
    return np.sqrt(hdot(x,x).real)
    
def uni_inte(h_or_hc, dt, c, opts={}):
    """
    unitary integration for Hermitian quantum system. 
    Compute time propagation of coefficient which satisfy below derivative equation
    .    idot{c} = H.c
    
    Inputs 
    ------
    h_or_hc : matrix or (vector->vector)
    .       Hamiltonian matrix H or function to calculate H.c

    InOut
    ------
    c       : vector
    """
    if("inte" not in opts):
        opts["inte"] = "diag"
    inte = opts["inte"]

    if("print_lvl" not in opts):
        opts["print_lvl"] = print_lvl = 0

    if(opts["print_lvl"]>0):
        print 
        print("math.uni_inte begin")
        print("inte:", opts["inte"])
        
    if(inte == "diag"):
        
        if(not isinstance(h_or_hc, np.ndarray)):
            raise RuntimeError("for method==diag, h_or_hc should be matrix")
        h = h_or_hc
        return uni_inte_diag(h, dt, c)
        
    elif(inte == "krylov"):
        strtype = str(type(h_or_hc))
        if(strtype == "<type 'function'>"):
            hc = h_or_hc            
        elif(isinstance(h_or_hc, np.ndarray)):
            hc = lambda c: dot(h_or_hc,c)
        else:
            raise RuntimeError("h_or_hc is invalid type")
        return uni_inte_krylov(hc, dt, opts["krylov_num"], c)

    elif(inte == "eig"):
        if(not isinstance(h_or_hc, tuple)):
            raise RuntimeError("method==eig need (e,u,uH) tuple")
        (e,u,uH) = h_or_hc
        return uni_inte_eig(e,u,uH,dt,c)
    
    else:
        raise RuntimeError("Unsupported inte.")
    
def uni_inte_eig(e, u, uH, dt, c):
    c = np.dot(uH, c)
    c = np.exp(-1.0j*e*dt) * c
    c = np.dot(u, c)
    return c
    
def uni_inte_diag(h, dt, c):
    (e,u) = np.linalg.eigh(h)
    uH = np.transpose(u.conj())
    c = uni_inte_eig(e,u,uH,dt,c)
    return c
    
def uni_inte_krylov(hc, dt, kn, c):

    # -- prepare memory --
    n = len(c)
    u = np.zeros((n,kn), dtype=complex)
    Hu = np.zeros((n,kn), dtype=complex)
    kh = np.zeros((kn,kn), dtype=complex)
    kc = np.zeros(kn, dtype=complex)
    ku = np.zeros((kn,kn), dtype=complex)
    kw = np.zeros(kn, dtype=float)

    # -- 1st process --
    u[:,0] = c
    u[:,0] /= norm(u[:,0])
    Hu[:,0] = hc(u[:,0])
    kh[0,0] = hdot(u[:,0], Hu[:,0])
    
    # -- 2nd process --
    u[:,1] = Hu[:,0] - kh[0,0]*u[:,0]
    u[:,1] /= norm(u[:,1])
    Hu[:,1] = hc(u[:,1])
    kh[1,1] = hdot(u[:,1], Hu[:,1])
    kh[0,1] = hdot(u[:,0], Hu[:,1])
    kh[1,0] = kh[0,1].conj()
    
    # -- proceeding process --
    for k in range(1,kn-1):
        u[:,k+1] = Hu[:,k] - kh[k-1,k]*u[:,k-1] -kh[k,k]*u[:,k]
        u[:,k+1] /= norm(u[:,k+1])
        Hu[:,k+1] = hc(u[:,k+1])
        kh[k+1,k+1] = hdot(u[:,k+1], Hu[:,k+1])
        kh[k,  k+1] = hdot(u[:,k],   Hu[:,k+1])
        kh[k+1,  k] = kh[k,  k+1].conj()
        
    # -- propagate --
    (kw, ku) = np.linalg.eigh(kh)
    ck = np.exp(-1.0j*kw*dt) * ku[0,:].conj()
    ck = dot(ku, ck)
    c = dot(u, ck)

    """
    for a in range(kn):
        for b in range(kn):
            print a,b,hdot(u[:,a], u[:,b]),hdot(u[:,a],Hu[:,b])
    """
    
    return c
    
