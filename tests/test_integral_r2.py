import math,unittest
from integral_r2 import reciprocal_points,weighted_box_integral_3d,integral_r2_contributions_3d

class TestIR2Boundary(unittest.TestCase):
    def test_infinity_mapping(self):
        q=reciprocal_points([(0,2,4)],(0,0,0))[0]
        self.assertTrue(math.isinf(q[0])); self.assertEqual(q[1],0.5)
    def test_degenerate_inf(self):
        self.assertEqual(weighted_box_integral_3d((math.inf,math.inf,0,1,0,1)),0.0)
    def test_boundary_contrib(self):
        c=integral_r2_contributions_3d([(0,2,2),(2,0,2)],(0,0,0))
        # normalized uniform simplex measure gives 2/9 for each contribution
        self.assertAlmostEqual(c[0],2/9,12); self.assertAlmostEqual(c[1],2/9,12)
if __name__=='__main__': unittest.main()
