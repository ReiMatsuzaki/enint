import numpy as np

from enint.molfunc import *

class GTO:
    def __init__(self, ex, c, R, n, katom=0, do_normalize=False):
        self.ncont = len(ex)
        self.ex = np.array(ex)
        self.c  = np.array(c)
        self.w = np.array(R)
        self.n = np.array(n)
        self.katom = 0
        if(do_normalize):
            smat = gtomat([self], op_s())
            self.c = self.c / np.sqrt(smat[0,0])
            
    def __str__(self):
        return """
        ==== GTO ====
        ex: {0}
        c:  {1}
        w:  {2}
        n:  {3}
        """.format(self.ex, self.c, self.w, self.n)

def prim_s(na, wa, za, nb, wb, zb):
    zp = za+zb
    wp = (za*wa + zb*wb)/zp
    d2 = sum([x*x for x in wa-wb])
    ep = np.exp(-za*zb/zp*d2)
    c = 1    
    for i in range(3):
        cd = coef_d(zp,wp[i],wa[i],wb[i],na[i],nb[i],0)
        c = c * cd
    return ep * (np.sqrt(np.pi/zp))**3 * c

def prim_r(i):
    def __func__(na, wa, za, nb, wb, zb):
        naa = [n for n in na]
        naa[i] = naa[i] + 1
        return (prim_s(naa,wa,za, nb,wb,zb)
            + wa[i] * prim_s(na,wa,za, nb,wb,zb))
    return __func__

def prim_dw(i):
    def __func__(na, wa, za, nb, wb, zb):
        #   Dw.(x-w)^n Exp[-z(x-w)^2]
        # = {-n(x-w)^{n-1} +2z(x-w)^{n+1}} Exp[-z(x-w)^2]
        nbp = [n for n in nb]
        nbp[i] = nbp[i]+1
        acc = 2*zb*prim_s(na,wa,za, nbp,wb,zb)
        if(nb[i]!=0):
            nbm = [n for n in nb]
            nbm[i] = nbm[i]-1
            acc = acc - nb[i]*prim_s(na,wa,za, nbm,wb,zb)
        return acc
    return __func__

def prim_na(wc):
    def __func__(na, wa, za, nb, wb, zb):
        zp = za+zb
        wp = (za*wa + zb*wb)/zp
        d = wp-wc
        d2p = dist2(d)
        d2  = dist2(wa-wb)
        ep  = np.exp(-za*zb*d2/zp)

        res = 0
        ns = np.zeros(3, dtype=int)
        for nx in range(na[0]+nb[0]+1):
            cx = coef_d(zp,wp[0],wa[0],wb[0],na[0],nb[0],nx)
            for ny in range(na[1]+nb[1]+1):
                cy = coef_d(zp,wp[1],wa[1],wb[1],na[1],nb[1],ny)
                for nz in range(na[2]+nb[2]+1):
                    ns[0]=nx; ns[1]=ny; ns[2]=nz; 
                    cz = coef_d(zp,wp[2],wa[2],wb[2],na[2],nb[2],nz)
                    cr = coef_R(zp,wp,wc,ns,0)
                    res += cx*cy*cz*cr
        return -2*np.pi*ep*res/zp
    return __func__

def prim_t(na, wa, za, nb, wb, zb):    
    
    s000 = prim_s(na, wa, za, nb, wb, zb)

    sp = [0,0,0]
    for i in range(3):
        npb = [n for n in nb]
        npb[i] = nb[i]+2
        sp[i] = prim_s(na, wa, za, npb, wb, zb)

    res = -2*zb*(2*nb[0] + 2*nb[1] + 2*nb[2] + 3) * s000
    res+=  4*zb*zb*(sp[0] + sp[1] + sp[2])

    for i in range(3):
        if(nb[i]>1):
            npb = [n for n in nb]
            npb[i] = nb[i]-2
            s = prim_s(na, wa, za, npb, wb, zb)
            res += nb[i]*(nb[i]-1) * s
    return -0.5*res

class OneOp():
    def __init__(self,prim,op_type="H"):
        self.op_type = op_type
        self.prim = prim

def op_s():
    return OneOp(prim_s)

def op_t():
    return OneOp(prim_t)

def op_r(i):
    return OneOp(prim_r)

def op_dw(i):
    return OneOp(prim_dw(i), "A")

def op_na(wc):
    return OneOp(prim_na(wc))

def gtoele(ga, op, gb):
    acc = 0.0
    for i in range(ga.ncont):
        for j in range(gb.ncont):
            acc += (ga.c[i]*gb.c[j]*
                    op.prim(ga.n, ga.w, ga.ex[i],
                            gb.n, gb.w, gb.ex[j]))
    return acc
    
def gtomat(gs, op):
    num = len(gs)
    mat = np.zeros((num,num))
    for a in range(num):
        for b in range(num):
            mat[a,b] = gtoele(gs[a], op, gs[b])
    return mat
            
def nshel2gto(j, normalize=False):
    ex = np.array(j["ex"])
    nshell = j["nshell"]
    ng = j["ng"]
    ian = j["ian"]
    c = np.array(j["c"])
    gs = []
    l_kmin_dict = {1:0, 2:1, 5:2, 11:3} # see inputa.src
    l_kmax_dict = {1:0, 4:1,10:2, 20:3} # see inputa.src
    for k in range(nshell):
        kstart = j["kstart"][k]
        katom =  j["katom"][k]
        ktype =  j["ktype"][k]
        kng   =  j["kng"][k]
        kex = [ex[kk-1] for kk in range(kstart,kstart+kng)]
        kks = range(kstart,kstart+kng)
        ls = range(l_kmin_dict[j["kmin"][k]],
                   l_kmax_dict[j["kmax"][k]]+1)
        kc = c[:,katom-1]
        if(0 in ls):
            kcs = [j["cs"][kk-1] for kk in kks]
            gs.append(GTO(kex, kcs, kc, [0,0,0], katom, normalize))
        if(1 in ls):
            kcp = [j["cp"][kk-1] for kk in kks]
            gs.append(GTO(kex, kcp, kc, [1,0,0], katom, normalize) )
            gs.append(GTO(kex, kcp, kc, [0,1,0], katom, normalize) )
            gs.append(GTO(kex, kcp, kc, [0,0,1], katom, normalize) )
        if(2 in ls):
            kcd = [j["cd"][kk-1] for kk in kks]
            gs.append(GTO(kex, kcd, kc, [2,0,0], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [0,2,0], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [0,0,2], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [1,1,0], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [1,0,1], katom, normalize))
            gs.append(GTO(kex, kcd, kc, [0,1,1], katom, normalize))
                      
    return gs


