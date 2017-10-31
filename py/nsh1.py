class Nshel_v1:
    def __init__(self):
        self.nshel = 0
        self.ng = 0
        
        self.ex = []        
        
        self.c  = [[], [], []]
        self.ian = []
        self.zan = []

        self.kstart = []
        self.katom  = []
        self.ktype  = []
        self.kng    = []
        self.kmin   = []
        self.kmax   = []
        self.kloc   = []

        self.ntable = np.array([
            [0,0,0],
            [0,0,0],
            [1,0,0],
            [0,1,0],
            [0,0,1],
            [2,0,0],
            [0,2,0],
            [0,0,2],
            [1,1,0],
            [1,0,1],
            [0,1,1]])
        
        self.ltable = np.sum(self.ntable, axis=1)

        num_n = len(self.ntable[:,0])
        self.coef = [[] for n in range(num_n)]
            
    def load(self, j):
        self.ex = np.array(j["ex"])
        self.nshel = j["nshell"]
        self.ng = j["ng"]
        self.ian = j["ian"]
        self.c  = np.array(j["c"])

        self.kstart = j["kstart"]
        self.katom  = j["katom"]
        self.ktype  = j["ktype"]
        self.kng    = j["kng"]
        self.kmin   = j["kmin"]
        self.kmax   = j["kmax"]
        self.kloc   = j["kloc"]

        cs = j["cs"]
        cp = j["cp"]
        cd = j["cd"]

        for k in range(len(self.ntable)):
            l = self.ltable[k]
            if(l==0):
                self.coef[k] = cs
            elif(l==1):
                self.coef[k] = cp
            elif(l==2):
                self.coef[k] = cd

    def add_atom(self, ci, anum, z):
        ia = len(self.c[0])
        for ir in range(3):
            self.c[ir].append(ci[ir])
        self.ian.append(anum)
        self.zan.append(z)
        return ia
    
    def add_shel(self, ltype, cs, ex, ia):
        self.ex += ex
        self.nshel += 1
        self.ng += len(ex)
            
        self.katom.append(ia)        
        self.ktype.append(0)
        self.kng.append(len(ex))

        
        if(isinstance(ltype, str)):
            if(ltype=="s"):
                self.kmin.append(1)
                self.kmax.append(1)
                l = 0                
            elif(ltype=="p"):
                self.kmin.append(2)
                self.kmax.append(4)
                l = 1
            elif(ltype=="d"):
                self.kmin.append(5)
                self.kmax.append(10)
                l = 2
            else:
                raise RuntimeError("not impl")
        elif(isinstance(ltype, list)):
            ins = self.get_ins(ltype)
            self.kmin.append(ins)
            self.kmax.append(ins)
            l = np.sum(ltype)
        else:
            raise RuntimeError("not impl")

        inlist = range(self.kmin[-1], self.kmax[-1]+1)
        num_n = len(self.ntable[:,0])
        c0 = [0 for n in range(len(cs))]
        for i in range(num_n):
            if(i in inlist):
                self.coef[i] += cs
            else:
                self.coef[i] += c0[:]

    def setup(self):

        self.kstart.append(1)
        self.kloc.append(1)
        for k in range(1,self.nshel):
            self.kstart.append(self.kstart[k-1] + self.kng[k-1])
            nloc = self.kmax[k-1] - self.kmin[k-1] + 1
            self.kloc.append(self.kloc[k-1] + nloc)

        self.c = np.array(self.c)

        self.coef = np.array(self.coef)
        coef_old = np.copy(self.coef)
        
        s = self.smat()
        idx = -1
        for ishel in range(self.nshel):
            ips = self.get_iprims(ishel)
            for k in range(self.kmin[ishel], self.kmax[ishel]+1):
                idx += 1
                for ip in ips:
                    self.coef[k,ip] = coef_old[k,ip]/np.sqrt(s[idx,idx])
                            
    def smat(self):
        mat = np.zeros((21,21))

        for ishel in range(self.nshel):
            ips = self.get_iprims(ishel)
            ns_i = self.ntable[self.kmin[ishel]:self.kmax[ishel]+1,:]
            max_ni = sum(ns_i[-1])
            iatom = self.katom[ishel]-1
            ci = self.c[:,iatom]
            
            
            for jshel in range(self.nshel):
                jps = self.get_iprims(jshel)
                ns_j = self.ntable[self.kmin[jshel]:self.kmax[jshel]+1,:]
                max_nj = sum(ns_j[-1])
                jatom = self.katom[jshel]-1
                cj = self.c[:,jatom]
                locj = self.kloc[jshel]-1
                
                d2 = sum([x*x for x in ci-cj])
                
                for ip in ips:
                    for jp in jps:
                        zi = self.ex[ip]
                        zj = self.ex[jp]                        
                        zp = zi+zj
                        wp = (zi*ci+zj*cj)/zp
                        ep = np.exp(-zi*zj/zp*d2)
                        ds = [ [ [ coef_d(zp,wp[ir],ci[ir],cj[ir],ni,nj,0)
                                   for nj in range(max_nj+1)]
                                 for ni in range(max_ni+1)]
                               for ir in range(3)]

                        i = self.kloc[ishel]-1-1
                        for ik in range(self.kmin[ishel],
                                       self.kmax[ishel]+1):
                            i += 1
                            j = self.kloc[jshel]-1-1
                            for jk in range(self.kmin[jshel],
                                       self.kmax[jshel]+1):
                                j += 1
                                acc = 1.0                                
                                for ir in range(3):
                                    ni = self.ntable[ik,ir]
                                    nj = self.ntable[jk,ir]
                                    acc *= ds[ir][ni][nj]
                                acc *= ep*(np.pi/zp)**(1.5)
                                li = self.ntable[ik]
                                lj = self.ntable[jk]
                                coef = self.coef[ik,ip]*self.coef[jk,jp]
                                mat[i,j]+=acc*coef
        return mat
                
    def get_iprims(self, ishel):
        start_i = self.kstart[ishel]
        ng_i = self.kng[ishel]
        ips = range(start_i-1,start_i+ng_i-1)
        return ips
        
    def get_ns(self, ishel):
        l_kmin_dict = {1:0, 2:1, 5:2, 11:3} # see inputa.src
        l_kmax_dict = {1:0, 4:1,10:2, 20:3} # see inputa.src
        ls = range(l_kmin_dict[self.kmin[ishel]],
                   l_kmax_dict[self.kmax[ishel]]+1)
        ns = []
        if(0 in ls):
            ns.append([0,0,0])
        if(1 in ls):
            ns.append([1,0,0])
            ns.append([0,1,0])
            ns.append([0,0,1])
        if(2 in ls):
            ns.append([2,0,0])
            ns.append([0,2,0])
            ns.append([0,0,2])
            ns.append([1,1,0])
            ns.append([1,0,1])
            ns.append([0,1,1])
        return (ns, np.max(ls))

    def get_ins(self, n):
        for ins in range(1, len(self.ntable[:,0])):
            if(n[0]==self.ntable[ins,0] and
               n[1]==self.ntable[ins,1] and
               n[2]==self.ntable[ins,2]):
                return ins
