from enint.nsh import *


print "mole_gammainc"
for (m,z) in [(0,1.0), (1,1.1), (2,3.3), (3,0.0)]:
    print m, z, mole_gammainc(m,z)

print ""
print "coef_R"
print coef_R_list(1.1, [0.0,0.1,0.2], [0.2,0.3,0.4], 1, 0)
