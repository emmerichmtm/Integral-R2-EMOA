import unittest
import numpy as np
from benchmarks import (
    ideal, pareto_reference, evaluate, n_var, ZDT3_IDEAL_F2,
    DTLZ7_A, DTLZ7_B, DTLZ7_C, DTLZ7_QMAX,
)


class TestExactBenchmarkIdeals(unittest.TestCase):
    def test_zero_ideals(self):
        for p in ('ZDT1','ZDT2'):
            np.testing.assert_allclose(ideal(p,2),(0.0,0.0),rtol=0,atol=0)
        for p in ('DTLZ1','DTLZ2'):
            np.testing.assert_allclose(ideal(p,3),(0.0,0.0,0.0),rtol=0,atol=0)

    def test_zdt3_exact_componentwise_ideal(self):
        self.assertAlmostEqual(ZDT3_IDEAL_F2,-0.7733690123266405,14)
        self.assertAlmostEqual(ideal('ZDT3',2)[1],ZDT3_IDEAL_F2,15)

    def test_dtlz7_exact_componentwise_ideal_3d(self):
        self.assertAlmostEqual(DTLZ7_A,0.25141183608891715,14)
        self.assertAlmostEqual(DTLZ7_B,0.6316265307000613,14)
        self.assertAlmostEqual(DTLZ7_C,0.8594008566447240,14)
        self.assertAlmostEqual(DTLZ7_QMAX,1.6929956344984225,14)
        z=ideal('DTLZ7',3)
        np.testing.assert_allclose(z,(0.0,0.0,2.614008731003155),rtol=0,atol=2e-14)

    def test_dtlz7_ideal_dominates_entire_sampled_pf(self):
        z=ideal('DTLZ7',3)
        pf=pareto_reference('DTLZ7',3,16000)
        self.assertTrue(np.all(pf >= z[None,:]-1e-13))
        # The f3 ideal is attained at the high-high corner of a PF patch.
        self.assertLessEqual(np.min(pf[:,2])-z[2],1e-12)

    def test_dtlz7_ideal_dominates_sampled_feasible_image(self):
        z=ideal('DTLZ7',3)
        rng=np.random.default_rng(1701)
        F=np.asarray([evaluate('DTLZ7',x,3)
                      for x in rng.random((2000,n_var('DTLZ7',3)))])
        self.assertTrue(np.all(F >= z[None,:]-1e-13))


if __name__=='__main__':
    unittest.main()
