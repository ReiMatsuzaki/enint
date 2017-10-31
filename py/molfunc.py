import numpy as np
from scipy.special import gamma, erf, gammainc

def prod(xs):
    acc = 1
    for x in xs:
       acc *= x
    return acc

def gtoint(n, z):
    if(n%2==1):
        return 0
    if(n==0):
        return np.sqrt(np.pi/z)
    return (n-1)/(2*z) * gtoint(n-2,z)

def mole_gammainc(m, z):
    """
    Compute incomplete gamma function defined below
    .    F_m(z) = Int_0^1 t^{2m} Exp[-zt^2]dt
    Impled with scipy.special.gammainc(a,x) defined below
    .    1/Gamma[a] Int_0^x t^{a-1} Exp[-t] dt
    Variable transformation 
    .    t->(t^2)/x
    .    dt -> (2t/x) dt
    .    gammainc(a,x) = x^{1-a} /Gamma[a]  Int_0^1 t^{2a-2} (2t/x) Exp[-(t^2)/x] dt
    .                  = 2 (x^(-a))/gamma[a] Int_0^1 t^{2a-1} Exp[-(t^2)/x] dt
    .   F_m(z) = (1/2) gamma[a] x^a 

    Special case
    .  z = 0
    .  F_m(0) = Int_0^1 t^{2m} dt 
    .         = (2m+1)^{-1} [t^{2m+1}]_0^1 
    .         = (2m+1)^{-1}
    """

    eps = 10.0**(-14.0)
    
    if(abs(z)<eps):
        return 1/(2*m+1.0)
    
    if(m==0):
        a = m+0.5
        res = gamma(a)/(2*z**a)*gammainc(a, z)
#        res = sqrt(pi/z)/2 * erf(sqrt(z))
    else:
        a = m+0.5
        res = gamma(a)/(2*z**a)*gammainc(a, z)
        
    if((not res<1) and (not res>-1)):
        raise RuntimeError("""gammain failed
        res: {2}
        m: {0}
        z: {1}
        """.format(m,z,res))
    return res
    
def coef_d(zp,wpk,wak,wbk,nak,nbk,nk):
    if(nak==0 and nbk ==0 and nk ==0):
        return 1.0
    if(nk<0 or nk>nak+nbk):
        return 0.0
    if(nak>0):
        return (1/(2*zp) * coef_d(zp,wpk,wak,wbk,nak-1,nbk,nk-1) +
                (wpk-wak)* coef_d(zp,wpk,wak,wbk,nak-1,nbk,nk)   +
                (nk+1.0) * coef_d(zp,wpk,wak,wbk,nak-1,nbk,nk+1)   )
    else:
        return (1/(2*zp) * coef_d(zp,wpk,wak,wbk,nak,nbk-1,nk-1) +
                (wpk-wbk)* coef_d(zp,wpk,wak,wbk,nak,nbk-1,nk)   +
                (nk+1.0) * coef_d(zp,wpk,wak,wbk,nak,nbk-1,nk+1)   )

def coef_R(zp,wpc,m,j):
    """
    Compute function R_{m,j}(zp,w,c). See T.Kuchitsu, J.Okuda and M.Tachikawa, Int.J.Q.Chem. 109, 540 (2008)

    Inputs
    -------
    zp: float
    wcp: [float]*3 
    .    vector Wp-Wc 
    m : [int]*3
    j : int

    Returns
    --------
    R_{m,j}(zp,w,c) : float
    """

    if(m[0]==0 and m[1]==0 and m[2]==0):
        try:
            return (-2*zp)**j * mole_gammainc(j,zp*dist2(wpc))
        except RuntimeError as e:
            raise RuntimeError("""
failed at mole_gammainc. 
zp = {0}
wp = {1}
c  = {2}
m  = {3}
j  = {4}
Below are errror message from mole_grammainc
{5}
            """.format(zp,wp,c,m,j,e.message))

    for i in range(3):
        im = np.zeros(3); im[i] = 1
        if(m[i]>0):
            res = wpc[i] * coef_R(zp,wpc,m-im,j+1)
            if(m[i]>1):
                res += (m[i]-1) * coef_R(zp,wpc,m-2*im,j+1)
            return res

    raise RuntimeError("one of m is negative integer")
        

def coef_d_list(zp,wp,wa,wb,mna,mnb):
    """
    ds = np.zeros((3,ma+1,mb+1,ma+mb+2))
    for i in range(3):
        for na in range(ma+1):
            for nb in range(mb+1):
                for n in range(ma+mb+2):
                    ds[i,na,nb,n] = coef_d(zp,wp[ir],wa[ir],wb[ir],nj,nk,n)
    """
    ds = [[[[coef_d(zp,wp[ir],wa[ir],wb[ir],nj,nk,n)
             for n in range(mna+mnb+2)]
            for nk in range(mnb+1)]
           for nj in range(mna+1)]
          for ir in range(3)]
    ds = np.array(ds)
    return ds

def coef_R_list(zp,wpc,maxn,n, method=1):
    if(method==0):
        rs = [[[coef_R(zp, wpc, [nx,ny,nz], n)
                for nz in range(maxn+1)]
               for ny in range(maxn+1)]
              for nx in range(maxn+1)]
        rs = np.array(rs)
        return rs
    elif(method==1):
        return coef_R_list_fast(zp,wpc,maxn)
    else:
        raise RuntimeError("not impl")

def coef_R_list_fast(zp,wpc,maxn):

    rmap = np.zeros((maxn+1,maxn+1,maxn+1,3*maxn+1))
    zwpc = zp*dist2(wpc)

    tmp = 1
    for j in range(3*maxn+1):
        rmap[0,0,0,j] = tmp * mole_gammainc(j,zwpc)
        tmp *= (-2*zp)

    for nnn in range(1,3*maxn+1):
        for nx in range(min(maxn, nnn)+1):
            for ny in range(min(maxn,nnn-nx)+1):
                for nz in range(min(maxn,nnn-nx-ny)+1):
                    for j in range(3*maxn-nnn+1):
                        if(nx>0):
                            tmp = wpc[0]*rmap[nx-1,ny,nz,j+1]
                            if(nx>1):
                                tmp += (nx-1)*rmap[nx-2,ny,nz,j+1]
                            rmap[nx,ny,nz,j] = tmp
                        elif(ny>0):
                            tmp = wpc[1]*rmap[nx,ny-1,nz,j+1]
                            if(ny>1):
                                tmp += (ny-1)*rmap[nx,ny-2,nz,j+1]
                            rmap[nx,ny,nz,j] = tmp
                        elif(nz>0):
                            tmp = wpc[2]*rmap[nx,ny,nz-1,j+1]
                            if(nz>1):
                                tmp += (nz-1)*rmap[nx,ny,nz-2,j+1]
                            rmap[nx,ny,nz,j] = tmp
    return rmap[:,:,:,0]
                        
def dist2(xs):
    acc = 0
    for x in xs:
        acc += x*x
    return acc
    
