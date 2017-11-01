import json
import numpy as np
from numpy import sqrt, exp, pi

from enint.molfunc import *
from enint.nsh0 import GTO

class PWShel:
    def __init__(self, ntypes, ex, r0, p0, g0):
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
        self.ex = ex
        self.cex = conjg(ex)
        self.coef = np.ones(num)
        self.ccoef = conjg(self.coef)

        self.r0 = np.array(r0)
        self.p0 = np.array(p0)
        self.g0 = g0
        self.j0 = 0

class PWNshel:
    def __init__(self):
        self.shelds = []

    def setup(self):
        jn = 0
        for shel in self.shels:
            shel.j0 = jn
            jn += shel.num
        s = self.smat()
        idx = -1
        for shel in self.shels:
            coef_old = np.copy(shel.coef)
            for jj in range(self.num):
                idx += 1
                shel.coef[jj,:]  = coef_old[jj,:]/np.sqrt(s[idx,idx])
                shel.ccoef[jj,:] = conjg(shel.coef[jj,:])
                
    def smat(self):
        n = self.num_basis()
        mat = np.zeros((n,n))

        for sj in self.shels:
            for sk in self.shels:
                rj = sj.r0
                rk = sk.r0
                pj = sj.p0
                pk = sk.p0
                gj = sj.g0
                gk = sk.g0
                d2 = sum([x*x for x in wj-wk])                

                for  jg in range(sj.ng):
                    for  kg in range(sk.ng):
                        zj = sj.cex[jg]
                        zk = sk.ex[kg]
                        zp = zj+zk
                        wp = (2*zj*rj + 2*zk*rk -1.0j*pj-1.0j*pk)/zp
                        ep = np.exp(-zj*rj*rj - zk*rk*rk
                                    +1.0j*rj*pj - 1.0j*rk*pk
                                    -1.0j*gj + 1.0j*gk
                                    +zp*wp**2)
                        cp = ep*(np.pi/zp)**(1.5)
                        ds = coef_d(zp,wp,wj,wk,sj.max_n,sk.max_n,0)

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

    def num_basis(self):
        return sum([shel.num for shel in self.shels])

    
                
