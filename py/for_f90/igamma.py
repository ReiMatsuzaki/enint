from enint.molfunc import igamma

for z in [0.1-0.2j, 0.3+0.1j, 21.2+15.5j, 23.0+15.5j]:
    calc = igamma(3, z)
    for m in range(4):
        print m, z, calc[m]
