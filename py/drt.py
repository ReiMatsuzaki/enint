import numpy as np

da = [0,0,1,1]
db = [0,1,-1,0]
dc = [1,0,1,0]
def make_abc_table(N,S,norbs):
    na = 1+N/2
    nb = 1+2*S+norbs
    nc = 1+norbs
    iabc_table = np.zeros((na,nb,nc), dtype=int)
    abc_table  = np.zeros((na*nb*nc,3), dtype=int)

    iabc = 0
    for a in range(na):
        for b in range(nb):
            for c in range(nc):
                iabc_table[a,b,c] = iabc
                abc_table[iabc, :] = [a,b,c]
                iabc += 1
    return (abc_table, iabc_table)

def correct_abc(a,b,c):
    return (a>=0 and b>=0 and c>=0)

class DRT_pre(object):
    def __init__(self, N, S, norbs):
        (self.abc_table, self.jabc_table) = make_abc_table(N,S,norbs)
        self.N = N
        self.S = S
        self.norbs = norbs
        self.k    = [[] for i in range(norbs+1)]
        self.jabc = [[] for i in range(norbs+1)]
        
        (a,b,c) = calc_abc(N,S,norbs)
        i = a+b+c        
        self.jabc[i] = [self.jabc_table[a,b,c]]
        self.add_k(i)

        for i in range(a+b+c-1, -1, -1):
            self.add_jabc(i)
            self.add_k(i)
        
    def show(self):
        print "   i  | jabc |  aj  bj  cj | k0j k1j k2j k3j"
        print "------+------+-------------+------------------"
        for i in range(self.norbs, -1, -1):
            for idx in range(len(self.jabc[i])):
                jabc = self.jabc[i][idx]
                [a,b,c] = self.abc_table[jabc,:]
                [k0,k1,k2,k3] = self.k[i][idx]
                if(idx==0):
                    line = " {0:3}  |".format(i)
                else:
                    line = "      |"
                line += " {0:3}  | {1:3} {2:3} {3:3} | {4:3} {5:3} {6:3} {7:3}" 
                print line.format(jabc, a, b, c, k0, k1, k2, k3)
        print "------+------+-------------+------------------"

    def add_k(self, i):
        for jabc in self.jabc[i]:
            k = [-1,-1,-1,-1]
            for d in range(4):
                aj = self.abc_table[jabc,0] - da[d]
                bj = self.abc_table[jabc,1] - db[d]
                cj = self.abc_table[jabc,2] - dc[d]
                if(correct_abc(aj,bj,cj)):
                    jabc_new = self.jabc_table[aj,bj,cj]
                    k[d] = jabc_new
            self.k[i].append(k)
                
    def add_jabc(self, i):

        for k in self.k[i+1]:
            for d in range(4):
                jabc = k[d]
                if(jabc!=-1 and jabc not in self.jabc[i]):
                    self.jabc[i].append(jabc)

def replace_if(v0, v1):
    def __func__(x):
        if(x==v0):
            return v1
        else:
            return x
    return __func__
    
class DRT(object):
    def __init__(self, N, S, norbs):
        self.N = N
        self.S = S
        self.norbs = norbs
        
        pre = DRT_pre(N,S,norbs)
        self.sort(pre)
        self.calc_weight()
        self.nwks = self.uxj[self.nj-1]

        self.uxj2lxj = np.zeros(self.nwks+1)
        for j in range(self.nj):
            self.uxj2lxj[self.uxj[j]] = self.lxj[j]

    def sort(self, pre):
        norbs = self.norbs
        nj = 0
        for i in range(norbs+1):
            nj += len(pre.jabc[i])

        jabc2j_table = -1+np.zeros(len(pre.abc_table[:,0]))
        j = 0
        for i in range(norbs, -1, -1):
            for jabc in pre.jabc[i]:
                jabc2j_table[jabc] = j
                j+=1

        self.j0 = np.zeros(norbs+1, dtype=int)
        self.j1 = np.zeros(norbs+1, dtype=int)
        self.abcj = np.zeros((nj, 3), dtype=int)
        self.kj   = np.zeros((nj, 4), dtype=int)
        j = 0
        self.j0[norbs] = 0
        for i in range(norbs, 0, -1):
            nDR = len(pre.k[i])
            self.j1[i]   = self.j0[i]+nDR-1
            self.j0[i-1] = self.j0[i]+nDR
            for (jabc, k) in zip(pre.jabc[i], pre.k[i]):
                self.abcj[j,:] = pre.abc_table[jabc,:]
                self.kj[j,:] = [jabc2j_table[kk] for kk in k]
                j += 1
                
        self.j0[0] = self.j1[1]+1
        self.j1[0] = self.j0[0]
        self.abcj[j,:] = 0
        self.kj[j,:] = -1
        self.nj = self.j1[0]+1

    def calc_weight(self):
        self.uxj = np.zeros(self.nj, dtype=int)        
        self.uxj[0] = 1
        for i in range(self.norbs, 0, -1):
            for j_i in range(self.j0[i],self.j1[i]+1):

                for d in range(4):
                    j_im = self.kj[j_i,d]
                    if(j_im != -1):
                        self.uxj[j_im] += self.uxj[j_i]

                    
        self.uyj = np.zeros((self.nj, 4), dtype=int)
        for i in range(self.norbs,0,-1):
            for j in range(self.j0[i], self.j1[i]+1):
                for d in range(4):
                    if(self.kj[j, d] != -1):
                        w = 0
                        for jj in range(self.j0[i], self.j1[i]+1):
                            for dd in range(4):
                                if(d>dd and self.kj[jj,dd] == self.kj[j,d]):
                                    w += self.uxj[jj]
                        self.uyj[j,d] = w
                    else:
                        self.uyj[j,d] = -1
        self.uyj[self.nj-1,:] = -1

        self.lxj = np.zeros(self.nj, dtype=int)
        self.lxj[-1] = 1
        for i in range(self.norbs):
            for j_ip in range(self.j0[i+1], self.j1[i+1]+1):
                for d in range(4):
                    j_i = self.kj[j_ip,d]
                    if(j_i != -1):
                        self.lxj[j_ip] += self.lxj[j_i]

        self.lyj = np.zeros((self.nj, 4), dtype=int)
        for i in range(1,self.norbs+1):
            for j in range(self.j0[i], self.j1[i]+1):
                w = 0
                for d in range(4):
                    if(self.kj[j,d] != -1):
                        self.lyj[j,d] = w
                        w += self.lxj[self.kj[j,d]]                        
                    else:
                        self.lyj[j,d] = -1
        self.lyj[self.nj-1,:] = -1
                        
    def build_1e_mat(self, m):
        pass
        """
        res = np.zeros((self.nwks, self.nwks))
        for p in range(self.norbs):
            for q in range(p):
                # q<p
        """
    def get_j(self, ds, i0):
        j = 0
        for i in range(self.norbs, i0-1, -1):
            j = self.kj[j,ds[i]]
        return j
        
    def next_uwk(self, ds_uwk):
        pass
    """
        for i in range(1, self.norbs):
            if(ds_uwk[i]!=3):
       """         
            
        
    def loop1(self, j_head, j_tail, ds_loop):
        if(not j_head<j_tail):
            raise RuntimeError("only j_head<j_tail")

        i_head = -1
        i_tail = -1
        for i in range(0, self.norbs+1):
            if(self.j0[i] <= j_head and j_head <= self.j1[i]):
                i_head = i
            if(self.j0[i] <= j_tail and j_tail <= self.j1[i]):
                i_tail = i

        """
        ds_last_uwk = np.zeros(self.norbs-i_head, dtype=int)
        for i in range(self.norbs, i_head-1, -1):
            for d in range(4):
                if(self.kj)
        """
            
    def show(self, l_or_u = "up"):
        if(l_or_u == "up"):
            print "   i  |   j  |  aj  bj  cj | k0j k1j k2j k3j |  uxj   uyj"
        elif(l_or_u == "low"):
            print "   i  |   j  |  aj  bj  cj | k0j k1j k2j k3j |  lxj   lyj"
        print "------+------+-------------+-----------------+----------------------"
        for i in range(self.norbs, -1, -1):
            j0 = self.j0[i]
            j1 = self.j1[i]
            for j in range(j0, j1+1):
                [a,b,c] = self.abcj[j,:]
                [k0,k1,k2,k3] = map(replace_if(-1,"  -"), self.kj[j,:])
                if(l_or_u == "up"):
                    [y0,y1,y2,y3] = map(replace_if(-1,"  -"), self.uyj[j,:])
                    xj = self.uxj[j]
                elif(l_or_u == "low"):
                    [y0,y1,y2,y3] = map(replace_if(-1,"  -"), self.lyj[j,:])
                    xj = self.lxj[j]
                if(j==j0):
                    line = " {0:3}  |".format(i)
                else:
                    line = "      |"
                line += " {0:3}  | {1:3} {2:3} {3:3} | {4:3} {5:3} {6:3} {7:3} | {8:3}  {9:3} {10:3} {11:3} {12:3}" 
                print line.format(j, a, b, c, k0, k1, k2, k3, xj,
                                  y0, y1, y2, y3)
        print "------+------+-------------+-----------------|---------------------"
                
class DRij(object):
    def __init__(self, a, b, c):
        self.abc = (a,b,c)
        self.nexts = []

    def show(self):
        (a,b,c) = self.abc
        print " {0:3}  | {1:3} {2:3} {3:3}".format(a+b+c,a,b,c)
        for dr in self.nexts:
            if(dr is not None):
                dr.show()

    def make_next(self):
        (a,b,c) = self.abc

        if(a+b+c==0):
            if(a==0 and b==0 and c==0):
                return True
            else:
                return False

        d0 = DRij(a,b,c-1)
        d1 = DRij(a,b-1,c)
        d2 = DRij(a-1,b+1,c-1)
        d3 = DRij(a-1,b,c)
        if(a==0 and b==0 and c==0):
            nexts = [None, None, None, None]
        elif(a==0 and b==0):
            nexts = [d0, None, None, None]
        elif(a==0):
            nexts = [d0, d1, None, None]
        elif(b==0):
            nexts = [d0, None, d2, d3]
        else:
            nexts = [d0, d1, d2, d3]

        res = False
        for n in nexts:
            if(n is None):
                res0 = False                
            else:
                res0 = n.make_next()
            if(res0):
                self.nexts.append(n)                
            else:
                self.nexts.append(None)
            res = res or res0
        return res

def calc_abc(N, S, n_orb):
    a = int(N/2.0-S)
    b = int(2*S)
    c = n_orb-a-b
    return (a, b, c)

def possible_ds(a, b, c):
    """
    di       |  0   1   2   3
    ---------+------------------
    delta ai |  0   0   1   1
    delta bi |  0   1  -1   0
    delta ci |  1   0   1   0
    """
    if(a==0 and b==0):
        return [0]
    elif(a==0):
        return [0,1]
    elif(b==0):
        return [0,2]
    else:
        return [0,1,2,3]

def calc_k(a, b, c, j):
    if(a==0 and b==0 and c==0):
        return None
    elif(a==0 and b==0):
        return [j+1,-1,-1,-1]
    elif(a==0):
        return [j+1,j+2,-1,-1]
    elif(b==0):
        return [j+1,-1,j+2,-1]
    else:
        return [j+1,j+2,j+3,j+4]
    
class DRT_old(object):
    def generate(self, i):
        """ from i+1 th (aj,bj,cj) and (k0j,k1j,k2j,k3j), compute i th DR.
        """
        
        if(i==0):
            self.j0[i] = self.j1[i+1]+1
            self.j1[i] = self.j1[i+1]+1
            j = self.j0[i]
            self.abc[j] = [0,0,0]
            self.k[j,:] = [-1,-1,-1,-1]
            return

        self.j0[i] = self.j1[i+1]+1
        j1 = self.j0[i]
        for j in range(self.j0[i+1], self.j1[i+1]+1):
            for d in range(4):
                if(self.k[j,d]!=-1):
                    if(j1<self.k[j,d]):
                        j1 = self.k[j,d]


        self.j1[i] = j1
        for j in range(self.j0[i+1], self.j1[i+1]+1):
            [a,b,c] = self.abc[j]
            for d in range(4):
                if(self.k[j,d]!=-1):
                    aa = a - da[d]
                    bb = b - db[d]
                    cc = c - dc[d]
                    kk = calc_k(aa,bb,cc,j1)
                    j1 += len([jk for jk in kk if jk!=-1])
                    self.abc[self.k[j,d],:] = [aa,bb,cc]
                    self.k[self.k[j,d],:] = kk

    def __init__(self, N, S, norbs, nalloc):
        self.N = N
        self.S = S
        self.norbs = norbs
        self.j0 = np.zeros(norbs+1, dtype=int)
        self.j1 = np.zeros(norbs+1, dtype=int)
        self.abc = np.zeros((nalloc, 3), dtype=int)
        self.k = np.zeros((nalloc, 4), dtype=int)

        self.j0[norbs] = 0
        self.j1[norbs] = 0
        (a,b,c) = calc_abc(N, S, norbs)
        self.abc[0,:] = [a,b,c]
        self.k[0,:] = calc_k(a,b,c,0)

        for i in range(norbs-1, -1, -1):
            self.generate(i)

    def __str__(self):
        
        res = "   i  |   j  |  aj  bj  cj | k0j k1j k2j k3j\n"
        res+= "------+------+-------------+------------------\n"
        for i in range(self.norbs, -1, -1):
            for j in range(self.j0[i], self.j1[i]+1):
                [a,b,c] = self.abc[j,:]
                [k0,k1,k2,k3] = self.k[j,:]
                if(j==self.j0[i]):
                    res += " {0:3}  |".format(i)
                else:
                    res += "      |"
                line = " {0:3}  | {1:3} {2:3} {3:3} | {4:3} {5:3} {6:3} {7:3}\n" 
                res += line.format(j, a, b, c, k0, k1, k2, k3)
        res+= "------+------+-------------+------------------\n"
        return res

    
