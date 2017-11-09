import os
join = os.path.join

import numpy as np
tr = np.transpose
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pylab as plt


from enint.nsh import *
from enint.math import ijv2mat
from enint.ciwfn import *

enint_root = os.path.abspath("../../")
qc_root = join(enint_root, "gms/h2/on_z/out")
nfrozen = 0

zs = np.linspace(-3.0,6.0,100)
rs = np.array([[0.0,0.0,z] for z in zs])

nsh = nshel_load(join(qc_root, "nshel.json"))
nsh.setup(True)
psi_ao = nsh.ao_at(rs)

cmo = ijv2mat(join(qc_root, "cmo.csv"))
rho0_mo = rho_mo(psi_ao, cmo)

aij = aij_load(join(qc_root, "aij.csv"), nfrozen=nfrozen)
cci = ijv2mat(join(qc_root, "cci.csv"))

for n in [0,1]:
    dm1 = aij.dm1(cci[:,n])
    rho = [ciwfn_op1(rho0_mo[:,:,ir], dm1) for ir in range(len(rs))]
    plt.plot(zs, rho, label=str(n))

plt.legend()
plt.savefig("rho.pdf")


