import numpy as np

da = [0,0,1,1]
db = [0,1,-1,0]
dc = [1,0,1,0]
nanval = -1
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

def calc_abc(N, S, n_orb):
    a = int(N/2.0-S)
    b = int(2*S)
    c = n_orb-a-b
    return (a, b, c)

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
            k = [nanval,nanval,nanval,nanval]
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
                if(jabc!=nanval and jabc not in self.jabc[i]):
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

        self.ydj = np.zeros((self.nj+1, 4), dtype=int)
        self.xj  = np.zeros(self.nj+1, dtype=int)
        self.u_ydj = np.zeros((self.nj+1, 4), dtype=int)
        self.u_xj  = np.zeros(self.nj+1, dtype=int)
        
        self.calc_weight()
        self.nwks = self.xj[self.nj-1]        

    def sort(self, pre):
        norbs = self.norbs
        nj = 0
        for i in range(norbs+1):
            nj += len(pre.jabc[i])
        self.nj = nj
        
        jabc2j_table = np.zeros(len(pre.abc_table[:,0]))
        j = 1
        for i in range(norbs, -1, -1):
            for jabc in pre.jabc[i]:
                jabc2j_table[jabc] = j
                j+=1
        
        self.j0 = np.zeros(norbs+1, dtype=int)
        self.j1 = np.zeros(norbs+1, dtype=int)
        self.abcj = np.zeros((nj+1, 3), dtype=int)
        self.kj   = np.zeros((nj+1, 4), dtype=int)
        j = 1
        self.j0[norbs] = j
        for i in range(norbs, 0, -1):
            nDR = len(pre.k[i])
            self.j1[i]   = self.j0[i]+nDR-1
            self.j0[i-1] = self.j0[i]+nDR
            for (jabc, k) in zip(pre.jabc[i], pre.k[i]):                
                self.abcj[j,:] = pre.abc_table[jabc,:]                
                for d in range(4):
                    if(k[d] != nanval):
                        self.kj[j,d] = jabc2j_table[k[d]]
                        
                    else:
                        self.kj[j,d] = nanval
                j += 1
                
        self.j0[0] = self.j1[1]+1
        self.j1[0] = self.j0[0]
        self.abcj[j,:] = 0
        self.kj[j,:] = nanval

        self.u_kj = np.zeros((nj+1, 4), dtype=int)
        self.u_kj[:,:] = nanval
        for i in range(norbs):
            for j_ip in range(self.j0[i+1], self.j1[i+1]+1):
                for d in range(4):
                    if(self.kj[j_ip, d] != nanval):
                        self.u_kj[self.kj[j_ip,d],d] = j_ip

    def calc_weight(self):
        self.xj[self.nj] = 1
        for i in range(1, self.norbs+1):
            for j in range(self.j0[i], self.j1[i]+1):
                dd = nanval
                for d in range(4):
                    if(self.kj[j,d] == nanval):
                        self.ydj[j,d] = nanval
                    else:
                        if(dd == nanval):
                            self.ydj[j,d] = 0
                            dd = d
                        else:
                            self.ydj[j,d] = self.ydj[j,dd] + self.xj[self.kj[j,dd]]
                            dd = d
                if(dd!=nanval):
                    self.xj[j] = self.ydj[j,dd] + self.xj[self.kj[j,dd]]
        
        self.u_xj[1] = 1
        for i in range(self.norbs, -1, -1):
            for j in range(self.j0[i], self.j1[i]+1):
                dd = nanval
                for d in range(4):
                    if(self.u_kj[j,d] == nanval):
                        self.u_ydj[j,d] = nanval
                    else:
                        if(dd == nanval):
                            self.u_ydj[j,d] = 0
                        else:
                            self.u_ydj[j,d] = (self.u_ydj[j,dd] +
                                               self.u_xj[self.u_kj[j,dd]])
                        dd = d
                if(dd!=nanval):
                    self.u_xj[j] = self.u_ydj[j,dd] + self.u_xj[self.u_kj[j,dd]]
                            
    def loop1(self, j_head, j_tail, ds_loop):
        if(not j_head<j_tail):
            raise RuntimeError("only j_head<j_tail")

        i_head = nanval
        i_tail = nanval
        for i in range(0, self.norbs+1):
            if(self.j0[i] <= j_head and j_head <= self.j1[i]):
                i_head = i
            if(self.j0[i] <= j_tail and j_tail <= self.j1[i]):
                i_tail = i

        last_uwk = np.zeros(self.norbs+1, dtype=int)
        last_uwk[:] = nanval
        j = j_head
        for i in range(i_head, self.norbs):
            for d in range(4):
                if(self.u_kj[j,d]!=nanval):
                    last_uwk[i+1] = d
                    j = self.u_kj[j,d]
                    break

        first_lwk = np.zeros(self.norbs+1, dtype=int)
        first_lwk[:] = nanval
        j = j_tail
        for i in range(i_tail, 0, -1):
            for d in range(4):
                if(self.kj[j,d]!=nanval):
                    first_lwk[i] = d
                    j = self.kj[j,d]
                    break

        ket_wk = (list(first_lwk[0:i_tail]) +
                  list(ds_loop[i_tail:i_head+1,1]) +
                  list(last_uwk[i_head+1:]))
        bra_wk = (list(first_lwk[0:i_tail]) +
                  list(ds_loop[i_tail:i_head+1,0]) +
                  list(last_uwk[i_head+1:]))

        j = 1
        ibra = 0
        for i in range(self.norbs,0,-1):            
            d = bra_wk[i]
            y = self.ydj[j,d]
            if(y==nanval):
                raise RuntimeError("not connect")
            ibra += self.ydj[j,d]
            j = self.kj[j,d]
        print ibra

        j = 1
        iket = 0
        for i in range(self.norbs,0,-1):
            d = ket_wk[i]
            iket += self.ydj[j,d]
            j = self.kj[j,d]
        print iket

        for nl in range(self.xj[j_tail]):
            "compute iket and ibra with upward"
                
        """
        ds_last_uwk = np.zeros(self.norbs-i_head, dtype=int)
        for i in range(self.norbs, i_head-1, -1):
            for d in range(4):
                if(self.kj)
        """
            
    def show(self, up_or_down = "down"):
        
        print "   i  |   j  |  aj  bj  cj | k0j k1j k2j k3j |  xj   ydj"
        print "------+------+-------------+-----------------+----------------------"
        for i in range(self.norbs, -1, -1):
            j0 = self.j0[i]
            j1 = self.j1[i]
            for j in range(j0, j1+1):
                [a,b,c] = self.abcj[j,:]                
                if(up_or_down == "down"):
                    [k0,k1,k2,k3] = map(replace_if(nanval,"  -"), self.kj[j,:])
                    [y0,y1,y2,y3] = map(replace_if(nanval,"  -"), self.ydj[j,:])
                    xj = self.xj[j]
                elif(up_or_down == "up"):
                    [k0,k1,k2,k3] = map(replace_if(nanval,"  -"), self.u_kj[j,:])
                    [y0,y1,y2,y3] = map(replace_if(nanval,"  -"), self.u_ydj[j,:])
                    xj = self.u_xj[j]
                if(j==j0):
                    line = " {0:3}  |".format(i)
                else:
                    line = "      |"
                line += " {0:3}  | {1:3} {2:3} {3:3} | {4:3} {5:3} {6:3} {7:3} | {8:3}  {9:3} {10:3} {11:3} {12:3}" 
                print line.format(j, a, b, c, k0, k1, k2, k3, xj,
                                  y0, y1, y2, y3)
        print "------+------+-------------+-----------------|---------------------"
                
