import unittest
import time
import os
import sys
join = os.path.join

import pandas as pd
import numpy as np
tr = np.transpose
dot = np.dot

from enint.nsh import *
from enint.ciwfn import *
from enint.math import ijv2mat

enint_root = os.path.abspath("../../")

class TestCIWfn(unittest.TestCase):
    def test_smat_h2(self):
        #qc_root = join(enint_root, "gms/lif_casscf66/naewdy2/out")
        qc_root = join(enint_root, "gms/h2/out")
        nel = 2
        nfrozen = 0
        
        nsh = nshel_load(join(qc_root, "nshel.json"))
        nsh.setup(True)
        sao = nsh.smat()

        cmo = ijv2mat(join(qc_root, "cmo.csv"))
        smo = ao2mo(sao, cmo)
        self.assertAlmostEqual(1.0, smo[0,0])
        self.assertAlmostEqual(0.0, smo[1,0])
        
        aij = aij_load(join(qc_root, "aij.csv"))
        smo /= nel
        scsf = aij.mo2csf(smo)
        self.assertAlmostEqual(1.0, scsf[0,0])
        self.assertAlmostEqual(0.0, scsf[1,0])
        
        cci = ijv2mat(join(qc_root, "cci.csv"))
        sci = aij.mo2ci(smo, cci)
        self.assertAlmostEqual(1.0, sci[0,0])
        self.assertAlmostEqual(0.0, sci[0,1])

        return

    def test_smat_LiH(self):
        qc_root = join(enint_root, "gms/lih/fz1/out")
        nel = 4
        
        nsh = nshel_load(join(qc_root, "nshel.json"))
        nsh.setup(True)
        sao = nsh.smat()

        cmo = ijv2mat(join(qc_root, "cmo.csv"))
        smo = ao2mo(sao, cmo)
        self.assertAlmostEqual(1.0, smo[0,0])
        self.assertAlmostEqual(0.0, smo[1,0])
        
        aij = aij_load(join(qc_root, "aij.csv"))
        smo /= nel
        scsf = aij.mo2csf(smo)
        self.assertAlmostEqual(1.0, scsf[0,0])
        self.assertAlmostEqual(0.0, scsf[1,0])
        
        cci = ijv2mat(join(qc_root, "cci.csv"))
        sci = aij.mo2ci(smo, cci)
        self.assertAlmostEqual(1.0, sci[0,0])
        self.assertAlmostEqual(0.0, sci[0,1])

        return
    
    def test_xdip_lif(self):
        print "aaaaaa"
        qc_root = join(enint_root, "gms/lif_casscf66/naewdy2/out")
        nel = 12
        
        nsh = nshel_load(join(qc_root, "nshel.json"))
        nsh.setup(True)
        sao = nsh.smat()
        xao = nsh.rmat(0)
        yao = nsh.rmat(1)
        zao = nsh.rmat(2)

        cmo = ijv2mat(join(qc_root, "cmo.csv"))
        smo = ao2mo(sao, cmo)
        xmo = ao2mo(xao, cmo)
        ymo = ao2mo(yao, cmo)
        zmo = ao2mo(zao, cmo)
        self.assertAlmostEqual(1.0, smo[0,0])
        self.assertAlmostEqual(0.0, smo[1,0])
        
        aij = aij_load(join(qc_root, "aij.csv"))
        smo /= nel
        scsf = aij.mo2csf(smo)
        xcsf = aij.mo2csf(xmo)

        self.assertAlmostEqual(1.0, scsf[0,0])
        self.assertAlmostEqual(0.0, scsf[1,0])
        
        cci = ijv2mat(join(qc_root, "cci.csv"))
        xci = -aij.mo2ci(xmo, cci) 
        yci = -aij.mo2ci(ymo, cci)
        zci = -aij.mo2ci(zmo, cci) + 3.0*9.0/3.0

        dm1 = aij.dm1(cci[:,0])
        zci1_dm = -ciwfn_op1(zmo, dm1) + 3.0*9.0/3.0
        
        dm1 = aij.dm1(cci[:,1])
        zci2_dm = -ciwfn_op1(zmo, dm1) + 3.0*9.0/3.0

        print "calc1(state1)", zci[0,0]
        print "calc2(state1)", zci1_dm
        print "ref(state1) ", 0.393456*6.671265
        print "calc1(state2)", zci[1,1]
        print "calc2(state2)", zci2_dm
        print "ref(state2)", 0.393456*(-3.538234)

        self.assertAlmostEqual(zci1_dm, zci[0,0])
        self.assertAlmostEqual(zci2_dm, zci[1,1])

        
        return
                
if __name__ == '__main__':
    unittest.main()
        
