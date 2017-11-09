import unittest
import time
import os
import sys
join = os.path.join

import pandas as pd
import numpy as np
tr = np.transpose
dot = np.dot

from enint.nsh0 import *
from enint.nsh import *
from enint.math import ijv2mat

enint_root = os.path.abspath("../../")

class TestCIWfn(unittest.TestCase):
    def test_first(self):
        qc_root = join(enint_root, "gms/lif_casscf66/naewdy2/out")
        
        nsh = nshel_load(join(qc_root, "nshel.json"))
        nsh.setup(True)
        xao = nsh.rmat(0)

        cmo = ijv2mat(join(qc_root, "cmo.csv"))
        xmo = ao2mo(xao, cmo)
        #xmo = dot(tr(cmo), dot(xao, cmo))

        cci = ijv2mat(join(qc_root, "cci.csv"))
        aij = aij_load(join(qc_root, "aij.csv"))
        xcsf = mo2csf(xmo, aij)
        xci = mo2ci(xmo, aij, cci)

        print xao
        print xmo
        
                
if __name__ == '__main__':
    unittest.main()
        
