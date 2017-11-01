import json
import numpy as np
from numpy import sqrt, exp, pi

from enint.molfunc import *
from enint.nsh0 import GTO

class Shel:
    def __init__(self, ntypes, ex, coef_l, w):

        if(isinstance(ntypes, str)):
            ntypes = [ntypes]

        self.ns = []
        if("s" in ntypes):
            self.ns.append([0,0,0])
        if("px" in ntypes):
            self.ns.append([1,0,0])            
        if("p" in ntypes):
            self.ns.append([1,0,0])
            self.ns.append([0,1,0])
            self.ns.append([0,0,1])
        if("dxx" in ntypes):
            self.ns.append([2,0,0])
        if("d" in ntypes):
            self.ns.append([2,0,0])
            self.ns.append([0,2,0])
            self.ns.append([0,0,2])
            self.ns.append([1,1,0])
            self.ns.append([1,0,1])
            self.ns.append([0,1,1])

        self.num = len(self.ns)

        self.ns = np.array(self.ns)
        self.max_n = max([sum(n) for n in self.ns])
        
        self.ex = ex
        self.ng = len(ex)

        self.coef = np.zeros((self.num, self.ng))
        for jn in range(self.num):
            nj = self.ns[jn]
            l = sum(nj)
            self.coef[jn,:] = coef_l[l]
        
        if((self.num, self.ng) != self.coef.shape):
            raise RuntimeError("size mismatch")

        self.w = np.array(w)

        self.j0 = 0
    
class Nucs:
    def __init__(self):
        self.ws = []
        self.anum = []
        self.zs = []
        self.num = 0
    def add_atom(self,w,anum,z):
        self.ws.append(w)
        self.anum.append(anum)
        self.zs.append(z)
        self.num += 1
        return self.num-1


class Nshel:
    def __init__(self, nucs):
        self.shels = []
        self.nucs = nucs

    def __str__(self):
        line = "jshel  w   n  ex  coef   \n"
        for jshel in range(len(self.shels)):
            shel = self.shels[jshel]
            
            line += " {0}  {1}\n".format(jshel, shel.w)

            for jn in range(shel.num):
                line += "    {0}\n".format(shel.ns[jn])
                for ig in range(shel.ng):
                    line += "     {0} {1} {2}\n".format("A", shel.ex[ig],
                                                        shel.coef[jn,ig])
        return line
            
    def add_shel(self,ntypes,ex,coef,ia):
        self.shels.append(Shel(ntypes,ex,coef,self.nucs.ws[ia]))
        
    def setup(self, normalize=True):
        jn = 0
        for shel in self.shels:
            shel.j0 = jn
            jn += shel.num

        if(normalize):
            s = self.smat()
            idx = -1
            for shel in self.shels:
                coef_old = np.copy(shel.coef)
                for jj in range(shel.num):
                    idx += 1
                    shel.coef[jj,:] = coef_old[jj,:]/np.sqrt(s[idx,idx])
    
    def smat(self):
        nn = sum([shel.num for shel in self.shels])
        mat = np.zeros((nn,nn))

        for sj in self.shels:
            for sk in self.shels:
                wj = sj.w
                wk = sk.w
                d2 = sum([x*x for x in wj-wk])                

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.ex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (zj*wj+zk*wk)/zp                        
                        ep = np.exp(-zj*zk/zp*d2)
                        cp = ep*(np.pi/zp)**(1.5)

                        ds = [ [ [ coef_d1(zp,wp[ir],wj[ir],wk[ir],nj,nk,0)
                                   for nk in range(sk.max_n+1)]
                                 for nj in range(sj.max_n+1)]
                               for ir in range(3)]

                        for jj in range(sj.num):
                            for kk in range(sk.num):
                                nj = sj.ns[jj]
                                nk = sk.ns[kk]
                                acc = prod([ds[ir][nj[ir]][nk[ir]]
                                            for ir in range(3)])
                                coef = cp * sj.coef[jj,jg] * sk.coef[kk,kg]
                                j = sj.j0 + jj
                                k = sk.j0 + kk
                                mat[j,k] += acc*coef
        return mat
    
    def tmat(self):
        nn = sum([shel.num for shel in self.shels])
        mat = np.zeros((nn,nn))

        for sj in self.shels:
            for sk in self.shels:
                wj = sj.w
                wk = sk.w
                d2 = sum([x*x for x in wj-wk])

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.ex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (zj*wj+zk*wk)/zp                        
                        ep = np.exp(-zj*zk/zp*d2)
                        cp = ep*(np.pi/zp)**(1.5)

                        ds = coef_d(zp,wp,wj,wk,sj.max_n,sk.max_n+2,0)

                        for jj in range(sj.num):
                            for kk in range(sk.num):
                                nj = sj.ns[jj]
                                nk = sk.ns[kk]
                                i3 = range(3)

                                acc = 0
                                s000 = prod([ds[i,nj[i],nk[i],0] for i in i3])
                                acc += -2*zk*(2*nk[0]+2*nk[1]+2*nk[2]+3)*s000
                                
                                for jr in i3:
                                    nkp = np.copy(nk[:])
                                    nkp[jr] += 2
                                    s = prod([ds[i,nj[i],nkp[i],0] for i in i3])
                                    acc += 4*zk*zk*s
                                    if(nk[jr]>1):
                                        nkm = np.copy(nk[:])
                                        nkm[jr] -= 2
                                        s = prod([ds[i,nj[i],nkm[i],0] for i in i3])
                                        acc += nk[jr]*(nk[jr]-1)*s
                                        
                                coef = cp * sj.coef[jj,jg] * sk.coef[kk,kg]
                                j = sj.j0 + jj
                                k = sk.j0 + kk
                                mat[j,k] += acc*coef
        mat *= -0.5
        return mat        
        
    def vmat(self, nucs=None):
        if(nucs==None):
            nucs = self.nucs
            
        nn = sum([shel.num for shel in self.shels])
        mat = np.zeros((nn,nn))

        for sj in self.shels:
            for sk in self.shels:
                wj = sj.w
                wk = sk.w
                d2 = sum([x*x for x in wj-wk])

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.ex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (zj*wj+zk*wk)/zp                        
                        ep = np.exp(-zj*zk/zp*d2)
                        cp = -2*np.pi*ep/zp
                        
                        ds = coef_d(zp,wp,wj,wk,sj.max_n,sk.max_n,sk.max_n+sj.max_n)
                        
                        for ic in range(nucs.num):
                            nmax = sj.max_n+sk.max_n
                            wpc = wp-nucs.ws[ic]
                            rs = coef_R(zp,wpc,nmax)
                              
                            for jj in range(sj.num):
                                for kk in range(sk.num):
                                    nj = sj.ns[jj]
                                    nk = sk.ns[kk]
                                    i3 = range(3)
                                    
                                    acc = 0
                                    for nx in range(nj[0]+nk[0]+1):
                                        for ny in range(nj[1]+nk[1]+1):
                                            for nz in range(nj[2]+nk[2]+1):
                                                acc += (ds[0,nj[0],nk[0],nx]*
                                                        ds[1,nj[1],nk[1],ny]*
                                                        ds[2,nj[2],nk[2],nz]*
                                                        rs[nx,ny,nz]*
                                                        nucs.zs[ic])
                                    coef = cp * sj.coef[jj,jg] * sk.coef[kk,kg]
                                    j = sj.j0 + jj
                                    k = sk.j0 + kk
                                    mat[j,k] += acc*coef
        return mat                

    def eri(self):
        pass
                
    def to_gtos(self):
        gtos = []
        for shel in self.shels:            
            for jn in range(shel.num):
                gtos.append(GTO(shel.ex,
                                shel.coef[jn,:],
                                shel.w,
                                shel.ns[jn,:]))
        return gtos

    def num_basis(self):
        nn = sum([shel.num for shel in self.shels])
        return nn
    
def nshel_load(j):
    
    zan = j["ian"]
    ian = j["ian"]
    w  = np.array(j["c"])
    natom = len(zan)
    nucs = Nucs()
    for ia in range(natom):
        nucs.add_atom(w[:,ia], ian[ia], zan[ia])

    nshel = Nshel(nucs)
    
    ex = np.array(j["ex"])
    num_shel = j["nshell"]
    ng = j["ng"]
    kstart = j["kstart"]
    katom  = j["katom"]
    ktype  = j["ktype"]
    kng    = j["kng"]
    kmin   = j["kmin"]
    kmax   = j["kmax"]
    kloc   = j["kloc"]

    cs = j["cs"]
    cp = j["cp"]
    cd = j["cd"]

    l_kmin = {1:0, 2:1, 5:2, 11:3} # see inputa.src
    l_kmax = {1:0, 4:1,10:2, 20:3} # see inputa.src

    for ishel in range(num_shel):
        ls = range(l_kmin[kmin[ishel]], l_kmax[kmax[ishel]]+1)
        spd = ["s", "p", "d"]
        ntypes = [spd[l] for l in ls]
        ig0 = kstart[ishel]-1
        ig1 = kstart[ishel] + kng[ishel]-1
        
        ex_i = ex[ig0:ig1]
        cs_i = cs[ig0:ig1]
        cp_i = cp[ig0:ig1]
        cd_i = cd[ig0:ig1]

        coef_l = [None,None,None]
        coef_l[0] = cs_i
        coef_l[1] = cp_i
        coef_l[2] = cd_i
        
        nshel.add_shel(ntypes, ex_i, coef_l, katom[ishel]-1)
        
    return nshel

        
