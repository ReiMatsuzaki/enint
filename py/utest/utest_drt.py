from enint.drt import *

nalloc = 100

"""
todo

rename and comment
i : level (0 to norbs)
j : Row number which identify the Distinct Row(DR)
kdj: Downward chaining indices or link index
ydj: Downward indexies, arc weight
xj : Downward indexies, node weight
u_kdj : Upward chaining indexies
u_ydj : Upward arc weight
u_xj  : Upward node weight

adaptation with fortran index

check definition of modify upward indexing. replation with level i.

"""

drt = DRT(4,0,4)
drt.show('down')
print
drt.show('up')

ds = np.zeros((4+1,2), dtype=int)
ds[:] = -1
ds[2,0] = 0
ds[1,0] = 3
ds[2,1] = 2
ds[1,1] = 1

drt.loop1(6, 13, ds)


