import numpy as np
from numpy import sqrt, pi, exp, cos, sin
from scipy.special import gamma, erf, gammainc, hyp1f1

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
def fact(n):
    acc = 1
    for i in range(n,0,-1):
        acc *= i
    return acc
    
def dfact(n):
    acc = 1
    for i in range(n,0,-2):
        acc *= i
    return acc

def comb(n, j):
    return fact(n)/(fact(j)*fact(n-j))
    
def igamma_f1(maxm, z):
    x = z.real
    y = z.imag
    eps = 10.0**(-10.0)

    if(x<eps or y<-eps):
        raise RuntimeError("Re[z] and Im[z] must be positive")

    z2 = x*x+y*y
    az = abs(z)
    c = sqrt(pi*(az+x)/(8*z2))
    s = sqrt(pi*(az-x)/(8*z2)) + eps

    nf_max = 50
    nf_0 = 10
    anR = np.zeros(nf_max+1)
    anI = np.zeros(nf_max+1)
    anR[0] = -x/(2*z2)*exp(-x)*cos(y) + y/(2*z2)*exp(-x)*sin(y)
    anI[0] = +y/(2*z2)*exp(-x)*cos(y) + x/(2*z2)*exp(-x)*sin(y)
    delta = 10.0**(-15)
    convq = False
    for n in range(1,nf_max):
        anR[n] = -(2*n-1) * (x/(2*z2)*anR[n-1] + y/(2*z2)*anI[n-1])
        anI[n] = +(2*n-1) * (y/(2*z2)*anR[n-1] - x/(2*z2)*anI[n-1])
        if(n>nf_0):
            if(c*delta > abs(anR[n]) and s*delta>abs(anI[n])):                
                convq = True
                break
            
    if(not convq):
        print 
        print c, delta, anR[nf_max-1]
        print c*delta > abs(anR[nf_max-1])
        print s, delta, anI[nf_max-1]
        print s*delta>abs(anI[nf_max-1])
        raise RuntimeError("""not converged.
z = {0}
c = {1}
s = {2}
delta = {3}
anR[NF_max] = {4}
        """.format(z, c, s, delta, anR[nf_max]))

    fmR = np.zeros(maxm+1)
    fmI = np.zeros(maxm+1)
    fmR[0] = +c + np.sum(anR)
    fmI[0] = -s + np.sum(anI)
    for m in range(1,maxm+1):
        fmR[m] = (2*m-1) * (x/(2*z2)*fmR[m-1] + y/(2*z2)*fmI[m-1])
        fmI[m] = (2*m-1) * (x/(2*z2)*fmI[m-1] - y/(2*z2)*fmR[m-1])

    bmR = np.zeros(maxm+1)
    bmI = np.zeros(maxm+1)
    bmR[0] = 0.0
    bmI[0] = 0.0
    for m in range(1,maxm+1):
        bmR[m] = anR[0] + (2*m-1)*(x/(2*z2)*bmR[m-1] + y/(2*z2)*bmI[m-1])
        bmI[m] = anI[0] + (2*m-1)*(x/(2*z2)*bmI[m-1] - y/(2*z2)*bmR[m-1])

    res = (fmR+bmR) + 1.0j*(fmI+bmI)
    return res

def igamma_f2(maxm, z):
    x = z.real
    y = z.imag
    eps = 10.0**(-10.0)

    if(x < -eps):
        raise RuntimeError("Re[z] must be positive")

    if(y < -eps):
        tmp = igamma_f2(maxm, x-1.0j*y)
        return tmp.conj()

    NR = 47
    bn = np.zeros(NR+1,dtype=complex)
    an = np.zeros(NR+1,dtype=complex)

    bn[0] = 1.0
    bn[1] = 1.0 + 0.5*z
    for n in range(2,NR+1):
        bn[n] = bn[n-1] + z*z/(4*(2*n-1)*(2*n-3))*bn[n-2]

    res = []
    for m in range(maxm+1):
        an[0] = 1.0
        t1 = float(2*m+1)/(2*m+3)
        t2 = float(2*m+1)/((2*m+3)*(2*m+5))
        t3 = float((2*m+1)*((2*m+1)**2+44))/(60*(2*m+3)*(2*m+5)*(2*m+7))
        an[1] = bn[1] - t1*z
        an[2] = bn[2] - t1*z - t2*z*z
        an[3] = bn[3] - t1*z - t2*z*z - t3*z*z*z
        for n in range(4,NR+1):
            f1 = float(2*n-2*m-5)/(2*(2*n-3)*(2*n+2*m+1))
            f2 = 1/float(4*(2*n-1)*(2*n-3))
            f3 = -f1/(4*(2*n-3)*(2*n-5))
            e = -f1
            an[n] = (1+f1*z)*an[n-1] + (e+f2*z)*z*an[n-2] + f3*z**3*an[n-3]
        res.append(1.0/(2*m+1) *an[NR]/bn[NR])
    return np.array(res)

def exp_igamma_g1(maxm, z):
    MAX_M = 100
    NF_MAX=50
    
    x = z.real
    y = z.imag
    eps = 10.0**(-10.0)
    delta = 10.0**(-15.0)
    

    if(maxm > MAX_M):
        raise RuntimeError("maxm must be lesser than 100")

    if(x < eps or y < -eps):
        raise RuntimeError("Re[z] and Im[z] must be positive")

    xxyy = x*x+y*y
    az = sqrt(xxyy)
    twp_xxyy = 2*xxyy

    c = sqrt((pi*(az+x)/(8*xxyy)))
    s = sqrt((pi*(az-x)/(8*xxyy)))
    an_R = 0.0
    an_I = 0.0
    anm1_R = 0.0
    anm1_I = 0.0
    conv = False
    snR = c*exp(-x)*sin(y) + s*exp(-x)*cos(y) + anm1_R
    snI = c*exp(-x)*cos(y) - s*exp(-x)*sin(y) + anm1_I
    for n in range(1,NF_MAX):
        an_R = (2*n-1) * (x/two_xxyy*anm1_R + y/two_xxyy*anm1_I)
        an_I = (2*n-1) * (x/two_xxyy*anm1_I - y/two_xxyy*anm1_R)
        sn_R += an_R
        sn_I += an_I
        anm1_R = an_R
        anm1_I = an_I
        if(delta>abs(an_R/sn_R) and delta>abs(an_I/sn_I)):
            conv = True
            break
    if(not convq):
        raise RuntimeError("not converged")

    g_R = np.zeros(MAX_M)
    g_I = np.zeros(MAX_M)
    g_R[0] = sn_R
    g_I[0] = sn_I
    for m in range(1,maxm+1):
        g_R[m] = -(2*m-1)*(x/two_xxyy*g_R[m-1] + y/two_xxyy*g_I[m-1])
        g_I[m] = +(2*m-1)*(y/two_xxyy*g_R[m-1] - x/two_xxyy*g_I[m-1])
    b_R = np.zeros(MAX_M)
    b_I = np.zeros(MAX_M)
    b_R[0] = 0.0
    b_I[0] = 0.0
    for m in range(1,maxm+1):
        b_R[m] = +x/two_xxyy - (2*m-1)*(x/two_xxyy*b_R[m-1] + y/two_xxyy*b_I[m-1])
        b_I[m] = -y/two_xxyy + (2*m-1)*(y/two_xxyy*b_R[m-1] - x/two_xxyy*b_I[m-1])

    return (g_R+b_R) + 1.0j*(g_I+b_I)

def exp_igamma_g2(maxm, z):
    MAX_N = 40
    x = z.real
    y = z.imag
    delta = 10.0**(-15.0)

    if(x<-delta):
        raise RuntimeError("Re[z] must be positive")

    res = []
    for m in range(maxm+1):
        conv = False
        
        bnm2_R=1.0
        bnm2_I=0.0
        bnm1_R=bnm2_R + 2.0/(2*m+5)*x
        bnm1_I=bnm2_I + 2.0/(2*m+5)*y

        anm2_R=1.0
        anm2_I=0.0
        anm1_R=anm2_R - 2*x/(2*m+3)
        anm1_I=anm2_I - 2*y/(2*m+3)

        for n in range(2,maxn):
            f1 = float(2*(2*m+1))/((4*n+2*m+1)*(4*n+2*m-3))
            f2 = (float(8*(n-1)*(2*n+2*m-1)) /
                  ((4*n+2*m-1)*(4*n+2*m-3)*(4*n+2*m-3)*(4*n+2*m-5)) )
            t1 = 1+f1*x
            t2 = f1*y
            t3 = (x*x-y*y)*f2
            t4 = 2*x*y*f2
            an_R = t1*anm1_R - t2*anm1_I + t3*anm2_R - t4*anm2_I
            an_I = t1*anm1_I + t2*anm1_R + t3*anm2_I + t4*anm2_R
            bn_R = t1*bnm1_R - t2*bnm1_I + t3*bnm2_R - t4*bnm2_I
            bn_I = t1*bnm1_I + t2*bnm1_R + t3*bnm2_I + t4*bnm2_R

            bn = 0.0
            for j in range(n+1):
                bn += (float(dfact(2*n+2*m+2*j+1))/dfact(4*n+2*m+1) *
                       comb(n, j) * (2*z)**(n-j))
            gn = (an_R+1.0j*an_I)/(bn_R+1.0j*bn_I)
            gnm1 = (anm1_R+1.0j*anm1_I)/(bnm1_R+1.0j*bnm1_I)
            if(abs((gn-gnm1)/gn) < delta):
                conv = True
            anm2_R=anm1_R; bnm2_R=bnm1_R; 
            anm2_I=anm1_I; bnm2_I=bnm1_I
            anm1_R=an_R;   bnm1_R=bn_R; 
            anm1_I=an_I;   bnm1_I=bn_I

        if(not conv):
            raise RuntimeError("not converged")

        res.append(1.0/(2*m+1) * (an_R+1.0j*an_I)/(bn_R+1.0j*an_I) )
    return np.array(res)
    
def igamma_py(maxm, z):
    res = [ 1.0/(2*m+1) * hyp1f1(m+0.5, m+1.5, -z)
            for m in range(0,maxm+1)]
    return np.array(res)

def igamma_real(maxm, z):
    
    eps = 10.0**(-14.0)
    res_list = []
    for m in range(maxm+1):
        if(abs(z)<eps):
            res_list.append(1/(2*m+1.0))
        else:
            a = m+0.5
            res = gamma(a)/(2*z**a)*gammainc(a, z)
            #        res = sqrt(pi/z)/2 * erf(sqrt(z))        
    
            if((not res<1) and (not res>-1)):
                raise RuntimeError("""gammain failed
                res: {2}
                m: {0}
                z: {1}
                """.format(m,z,res))
            res_list.append(res)
    return np.array(res_list)
    
def igamma(maxm, z):
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

    x = z.real
    y = z.imag
    if(isinstance(z, float)):
        return in_gamma_real(maxm, z)
    elif(isinstance(z, complex)):
        eps = 10.0**(-15)
        if(x > -eps and y > -eps):
            if(x < 21 and x+y < 37):
                return igamma_f2(maxm, z)
            else:
                return igamma_f1(maxm, z)
        else:
            res = igamma(maxm, z.conjugate())
            return res.conjugate()
            
def coef_d1(zp,wpk,wak,wbk,nak,nbk,nk):
    if(nak==0 and nbk ==0 and nk ==0):
        return 1.0
    if(nk<0 or nk>nak+nbk):
        return 0.0
    if(nak>0):
        return (1/(2*zp) * coef_d1(zp,wpk,wak,wbk,nak-1,nbk,nk-1) +
                (wpk-wak)* coef_d1(zp,wpk,wak,wbk,nak-1,nbk,nk)   +
                (nk+1.0) * coef_d1(zp,wpk,wak,wbk,nak-1,nbk,nk+1)   )
    else:
        return (1/(2*zp) * coef_d1(zp,wpk,wak,wbk,nak,nbk-1,nk-1) +
                (wpk-wbk)* coef_d1(zp,wpk,wak,wbk,nak,nbk-1,nk)   +
                (nk+1.0) * coef_d1(zp,wpk,wak,wbk,nak,nbk-1,nk+1)   )

def coef_R1(zp,wpc,m,j):
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
            res = wpc[i] * coef_R1(zp,wpc,m-im,j+1)
            if(m[i]>1):
                res += (m[i]-1) * coef_R1(zp,wpc,m-2*im,j+1)
            return res

    raise RuntimeError("one of m is negative integer")
        

def coef_d(zp,wp,wa,wb,mna,mnb,maxn):
    ds = [[[[coef_d1(zp,wp[ir],wa[ir],wb[ir],nj,nk,n)
             for n in range(maxn+1)]
            for nk in range(mnb+1)]
           for nj in range(mna+1)]
          for ir in range(3)]
    ds = np.array(ds)
    return ds

def coef_R(zp,wpc,maxn, method=1):
    if(method==0):
        rs = [[[coef_R1(zp, wpc, [nx,ny,nz], 0)
                for nz in range(maxn+1)]
               for ny in range(maxn+1)]
              for nx in range(maxn+1)]
        rs = np.array(rs)
        return rs
    elif(method==1):
        return coef_R_fast(zp,wpc,maxn)
    else:
        raise RuntimeError("not impl")

def coef_R_fast(zp,wpc,maxn):

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
    
