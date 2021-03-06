import unittest
import time
import os
import sys
join = os.path.join

import pandas as pd

from enint.nsh0 import *
from enint.nsh import *
from enint.math import ijv2mat

class TestCase(unittest.TestCase):
    def assertMatProp(self,mattype,a,msg=""):
        (n,m) = a.shape
        mattype_list = [mattype]
        if(mattype=="overlap"):
            mattype_list.append("hermite")

        if "overlap" in mattype_list:
            self.assertEqual(n,m,msg=msg+"\noverlap must be squared")
            for i in range(n):
                self.assertAlmostEqual(1.0, a[i,i],
                                       msg="""
{2}
diagonal element need to unity
i      : {0}
a[i,i] : {1}
""".format(i, a[i,i], msg))

        if "hermite" in mattype_list:
            for i in range(n):
                for j in range(i):
                    self.assertAlmostEqual(a[i,j], a[j,i], msg="""
{4}
hermiticity is broken
i  : {0}
j  : {1}
a[i,j] : {2}
a[j,i] : {3}
""".format(i,j,a[i,j],a[j,i],msg))
        
    def assertMatEqual(self,a,b,prec=7,msg=""):
        self.assertEqual(a.shape,b.shape, msg="""{2}                         
shape of a and b are different
a: {0}
b: {1}
        """.format(a.shape, b.shape, msg))
        
        (n,m) = a.shape
        for i in range(n):
            for j in range(m):
                self.assertAlmostEqual(a[i,j], b[i,j], prec, msg="""{4}
Matrix a and b is largely different at (i,j).
(i,j) = ({0},{1})
a[i,j] = {2}
b[i,j] = {3}
""".format(i,j,a[i,j],b[i,j],msg))

    
class TestNsh(TestCase):
    def test_dz(self):
        
        g1 = GTO([1.5, 0.3], [1.0, 0.8], [0.0,0.0,0.0], [2,0,0], True)
    
        ex = [1.3, 0.5]
        cs = [0.8, 0.4]
        r  = np.array([0.0, 0.0, 1.1])
        dx = 0.001
        dr = np.array([0.0, 0.0, dx])
        n  = [0,0,2]
        g0 = GTO(ex, cs, r,    n, True)
        gp = GTO(ex, cs, r+dr, n, True)
        gm = GTO(ex, cs, r-dr, n, True)

        op = op_s()
        ref = (gtoele(g1,op_s(),gp) - gtoele(g1,op_s(),gm))/(2*dx)
        calc = gtoele(g1, op_dw(2), g0)
        self.assertAlmostEqual(ref, calc, 5)
        
    def test_dz2(self):

        ex = [1.3, 0.5]
        cs = [0.8, 0.4]
        r  = np.array([0.0, 0.0, 1.1])
        dx = 0.001
        dr = np.array([0.0, 0.0, dx])
        
        nucs = Nucs()
        nucs.add_atom([0.0,0.0,0.0],  1, 1.0)
        nucs.add_atom(r, 2, 1.0)
        nucs.add_atom(r+dr, 2, 1.0)
        nucs.add_atom(r-dr, 2, 1.0)
        nsh = Nshel(nucs)
        nsh.add_shel("dxx", [1.5, 0.3], {2:[1.0, 0.8]}, 0)
        nsh.add_shel("dzz", ex, {2:cs}, 1)
        nsh.add_shel("dzz", ex, {2:cs}, 2)
        nsh.add_shel("dzz", ex, {2:cs}, 3)
        
        nsh.setup(True)
        ref   = (nsh.smat()[0,2]-nsh.smat()[0,3])/(2*dx)
        calc  = nsh.dwmat(2)[0,1]        
        self.assertAlmostEqual(ref, calc, 7)
        
    def test_igamma(self):
        from scipy import integrate
        import numpy as np

        t = np.linspace(0,1.0,220)
        for (m,z) in [(0,1.0), (0,2.0), (1,1.1), (2,1.2), (4,1.1), (1,0.0)]:
            y = t**(2*m) * np.exp(-z*t*t)
            ref = integrate.simps(y, t)
            calc = igamma(m, z)[m]
            self.assertAlmostEqual(ref, calc,
                                   msg="""
ref = {0}
calc = {1}
|ref-calc| = {4}
(m,z) = ({2},{3})""".format(ref,calc,m,z,abs(ref-calc)))
                    
    def test_coef_R(self):
        wp = np.array([0.0, 0.1, 0.2])
        wc = np.array([0.2, 0.3, 0.4])
        wpc = wp-wc
        
        rs = coef_R(1.1, wpc, 1, method=0)
        self.assertAlmostEqual(0.95768901, rs[0,0,0])
        self.assertAlmostEqual(0.0352501,  rs[1,1,0])
        self.assertAlmostEqual(0.0352501,  rs[0,1,1])

        rs = coef_R(1.1, wpc, 1, 1)
        self.assertAlmostEqual(0.95768901, rs[0,0,0])
        self.assertAlmostEqual(0.0352501,  rs[1,1,0])
        self.assertAlmostEqual(0.0352501,  rs[0,1,1])
        
    def test_nshel(self):

        nucs = Nucs()
        
        ia1 = nucs.add_atom([0.0,0.0,0.0], 1, 1.0)
        ia2 = nucs.add_atom([1.0,0.0,0.0], 2, 0.3)        
        
        nshel = Nshel(nucs)
        nshel.add_shel("s", [1.1],   {0:[1.0]}, ia1)
        nshel.add_shel("p", [1.3],   {1:[1.0]}, ia2)
        nshel.add_shel("dxx", [1.2], {2: [1.0]}, ia2)

        self.assertEqual(ia1, nshel.shels[0].ia)
        self.assertEqual(ia2, nshel.shels[1].ia)
        self.assertEqual(ia2, nshel.shels[2].ia)
        
        nshel.setup(True)
        gs = nshel.to_gtos()

        ref  = np.array([ia1,ia2,ia2,ia2,ia2])
        calc = nshel.ia_vec()
        self.assertEqual(ref[0], calc[0])
        self.assertEqual(ref[1], calc[1])
        self.assertEqual(ref[2], calc[2])
        self.assertEqual(ref[3], calc[3])
        self.assertEqual(ref[4], calc[4])
        
        calc = nshel.smat()
        ref = gtomat(gs, op_s())
        self.assertMatProp("overlap", calc)
        self.assertMatEqual(ref, calc)

        calc = nshel.tmat()
        ref = gtomat(gs, op_t())
        self.assertMatProp("hermite", calc)
        self.assertMatEqual(ref, calc)

        calc = nshel.vmat()
        ref = gtomat(gs, op_na(nucs.ws[0])) + 0.3*gtomat(gs, op_na(nucs.ws[1]))
        self.assertMatProp("hermite", calc)
        self.assertMatEqual(ref, calc, msg="test_nshel.Check Nuclear Attraction")
        
    def test_nshel_rmat(self):
        nucs = Nucs()
        ia1 = nucs.add_atom([0.1,0.2,0.3], 1, 1.0)
        ia2 = nucs.add_atom([0.0,0.0,0.0], 2, 2.0)
        nshel = Nshel(nucs)
        nshel.add_shel("s", [1.2], {0:[1.4]}, ia1)
        nshel.add_shel("s", [1.1], {0:[1.3]}, ia2)
        nshel.add_shel("p", [1.1], {1:[1.3]}, ia2)
        nshel.setup(False)

        xmat = nshel.rmat(0)
        ymat = nshel.rmat(1)
        zmat = nshel.rmat(2)
        smat = nshel.smat()

        self.assertMatProp("hermite", xmat)
        self.assertMatProp("hermite", ymat)
        self.assertMatProp("hermite", zmat)

        self.assertAlmostEqual(0.0, xmat[1,1])
        self.assertAlmostEqual(0.0, ymat[1,1])
        self.assertAlmostEqual(0.0, zmat[1,1])

        self.assertAlmostEqual(smat[0,2], xmat[0,1])
        self.assertAlmostEqual(smat[0,3], ymat[0,1])
        self.assertAlmostEqual(smat[0,4], zmat[0,1])
        
    def test_nshel_h2(self):
        out = "../../gms/h2/out"
        with open(os.path.join(out, "nshel.json")) as f:
            j = json.load(f)
            nshel = nshel_load(j)
            nshel.setup()
            
        calc = nshel.smat()
        df = pd.read_csv(os.path.join(out, "s.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_nshel_h2")

        calc = nshel.tmat()
        df = pd.read_csv(os.path.join(out, "t.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_nshel_h2. T matrix")

        calc = nshel.tmat() + nshel.vmat()
        df = pd.read_csv(os.path.join(out, "h.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_nshel_h2. H core matrix")        
        
    def test_nshel_gms(self):
        out = "../../gms/hcp/out"

        with open(join(out, "nshel.json")) as f:
            j = json.load(f)
            nshel = nshel_load(j)
            nshel.setup(True)


        calc = nshel.smat()
#        self.assertMatProp("overlap", calc)
        df = pd.read_csv(join(out, "s.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc)

        calc = nshel.tmat()
#        self.assertMatProp("hermite", calc)
        df = pd.read_csv(join(out, "t.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc)

        calc = nshel.tmat() + nshel.vmat()
#        self.assertMatProp("hermite", calc)
        df = pd.read_csv(join(out, "h.csv"))
        ref  = ijv2mat(df)
        self.assertMatEqual(ref, calc, msg="test_neshl_gms. H core")        

    def test_ao_at(self):
        nucs = Nucs()
        d = 1.0
        ia1 = nucs.add_atom([0.0,0.0,0.0], 1, 1.0)
        ia2 = nucs.add_atom([d,0.0,0.0], 2, 0.3)
        
        nshel = Nshel(nucs)
        nshel.add_shel("s",   [1.1],   {0:[1.0]}, ia1)
        nshel.add_shel("p",   [1.3],   {1:[1.0]}, ia2)
        nshel.add_shel("dxx", [1.2], {2: [1.0]}, ia2)
        
        nshel.setup(True)
        
        rs = np.array([[0.0,0.0,0.0], [1.3,0.1,-0.1]])
        ao_rs = nshel.ao_at(rs, method=0)
        self.assertAlmostEqual(nshel.shels[0].coef[0,0], ao_rs[0,0])
        self.assertAlmostEqual(-nshel.shels[1].coef[0,0]*d*exp(-1.3*d*d), ao_rs[1,0])
        self.assertAlmostEqual(0.0, ao_rs[2,0])
        self.assertAlmostEqual(0.0, ao_rs[3,0])
        self.assertAlmostEqual(nshel.shels[2].coef[0,0]*d**2*exp(-1.2*d*d),
                               ao_rs[4,0])

        ao_rs1 = nshel.ao_at(rs, method=1)
        self.assertMatEqual(ao_rs, ao_rs1)

    def _test_ao_at_d(self):
        nucs = Nucs()
        ia0 = nucs.add_atom([0.0,0.1,0.2], 1, 0.3)
        ia1 = nucs.add_atom([0.1,0.2,0.1], 2, 0.4)
        ia2 = nucs.add_atom([0.0,0.4,0.1], 3, 0.5)
        
        nshel = Nshel(nucs)
        nshel.add_shel("s", [1.1], {0:[1.0]}, ia0)
        nshel.add_shel("p", [1.3], {1:[1.0]}, ia1)
        nshel.add_shel("d", [1.2], {2:[1.0]}, ia2)
        nshel.setup(True)

        r1 = np.array([0.1,0.2,0.3])
        dx = 0.01
        dr = dx * np.identity(3)
        drx = dr[:,0]
        dry = dr[:,1]
        drz = dr[:,2]
        rs = np.array([r1, r1+drx, r1-drx, r1+dry, r1-dry, r1+drz, r1-drz])

        phi0 = nshel.ao_at(rs, method=0)
        phi1 = nshel.ao_at(rs, method=1)
        self.assertMatEqual(phi0, phi1)

        phi_x = nshel.ao_at(rs, n_dr=[1,0,0], method=1)
        phi_y = nshel.ao_at(rs, n_dr=[0,1,0], method=1)
        phi_z = nshel.ao_at(rs, n_dr=[0,0,1], method=1)

        for mu in range(nshel.num_basis()):
            self.assertAlmostEqual(phi_x[mu,0], (phi0[mu,1]-phi0[mu,2])/(2*dx))
            self.assertAlmostEqual(phi_y[mu,0], (phi0[mu,3]-phi0[mu,4])/(2*dx))
            self.assertAlmostEqual(phi_z[mu,0], (phi0[mu,5]-phi0[mu,6])/(2*dx))
        
        
class TestCnshel(TestCase):
    def test_first(self):
        print coef_d1(1.0-0.1j,0.2,0.3,0.1,1,0,0)
        
if __name__ == '__main__':
    unittest.main()
        
