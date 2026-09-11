import unittest
from performance_indicators import ir2_value
from integral_r2 import integral_r2_contributions

class TestPerformanceIR2(unittest.TestCase):
    def test_singleton_2d(self):
        self.assertAlmostEqual(ir2_value([(1.0,1.0)],(0.0,0.0)),0.75,12)
    def test_singleton_3d(self):
        self.assertAlmostEqual(ir2_value([(1.0,1.0,1.0)],(0.0,0.0,0.0)),11.0/18.0,12)
    def test_contribution_matches_indicator_difference_2d(self):
        P=[(0.4,1.4),(0.7,0.9),(1.3,0.5)]
        z=(0.0,0.0)
        c=integral_r2_contributions(P,z)
        base=ir2_value(P,z)
        for i in range(len(P)):
            Q=P[:i]+P[i+1:]
            self.assertAlmostEqual(c[i],ir2_value(Q,z)-base,11)
    def test_contribution_matches_indicator_difference_3d(self):
        P=[(0.5,1.2,1.4),(0.8,0.8,1.1),(1.2,0.6,0.8),(1.5,1.0,0.5)]
        z=(0.0,0.0,0.0)
        c=integral_r2_contributions(P,z)
        base=ir2_value(P,z)
        for i in range(len(P)):
            Q=P[:i]+P[i+1:]
            self.assertAlmostEqual(c[i],ir2_value(Q,z)-base,10)

if __name__=='__main__': unittest.main()
