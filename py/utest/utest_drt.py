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

more simple code!
see Shavitt .ppt
"""

drt = DRT(4,0,4)
drt.show('up')
print
drt.show('low')

ds = np.zeros(4+1)
ds[4] = 3
ds[3] = 0
#print drt.get_j(ds, 2)

