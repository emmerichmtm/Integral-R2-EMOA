import math, unittest
import numpy as np
from benchmarks import ideal, ZDT3_IDEAL_F2, DTLZ7_QMAX

class TestExactBenchmarkIdeals(unittest.TestCase):
    def test_zero_ideals(self):
        for p in ('ZDT1','ZDT2','DTLZ1','DTLZ2'):
            np.testing.assert_allclose(ideal(p,2),(0.0,0.0),rtol=0,atol=0)

    def test_zdt3_exact_componentwise_ideal(self):
        self.assertAlmostEqual(ZDT3_IDEAL_F2,-0.7733690123266405,14)
        self.assertAlmostEqual(ideal('ZDT3',2)[1],ZDT3_IDEAL_F2,15)

    def test_dtlz7_exact_componentwise_ideal(self):
        self.assertAlmostEqual(DTLZ7_QMAX,1.6929956344984225,14)
        self.assertAlmostEqual(ideal('DTLZ7',2)[1],2.3070043655015775,14)

if __name__=='__main__': unittest.main()
