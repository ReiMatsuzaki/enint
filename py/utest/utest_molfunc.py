import unittest

from enint.molfunc import *

class TestMolfunc(unittest.TestCase):
    def test_fact(self):
        self.assertEqual(1, fact(1))
        self.assertEqual(2, fact(2))
        self.assertEqual(6, fact(3))
        self.assertEqual(24, fact(4))
        self.assertEqual(120, fact(5))

    def test_dfact(self):
        self.assertEqual(1, dfact(1))
        self.assertEqual(2, dfact(2))
        self.assertEqual(3, dfact(3))
        self.assertEqual(8, dfact(4))
        self.assertEqual(15, dfact(5))

    def test_comb(self):
        self.assertEqual(3, comb(3,1))
        self.assertEqual(3, comb(3,1))
        self.assertEqual(10, comb(5,2))
        self.assertEqual(10, comb(5,3))

    def test_igamma(self):
        ref = 0.96392503253589557+0.06262986651323893j
        
        calc = igamma_py(0, 0.1-0.2j)        
	self.assertAlmostEqual(ref, calc[0])

        calc = igamma_f2(0, 0.1-0.2j)
        self.assertAlmostEqual(ref, calc[0])

        for z in [0.1-0.2j, 0.3+0.1j, 21.2+15.5j, 23.0+15.5j]:
            res = igamma_py(3, z)
            calc = igamma(3, z)
            for i in range(3+1):
                self.assertAlmostEqual(res[i], calc[i])
        
if __name__ == '__main__':
    unittest.main()
        

