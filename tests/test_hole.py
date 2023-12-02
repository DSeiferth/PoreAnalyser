import unittest
#import PoreFinding_pdb as pf 
import sys
#sys.path.append('../')
import hole_analysis as hole_analysis
import numpy as np

class HoleTest(unittest.TestCase):
    """
    Tests the hole_analysis function.
    """
    def test_hole_analysis(self):
        """
        Tests hole_analysis function.
        """
        path = 'pdb_models/'
        name = '7tvi.pdb'
        #midpoints2, means2 = hole_analysis.hole_analysis(name, path, end_radius=20, sel='protein')

        #self.assertEqual(len(midpoints2), len(means2))
        #self.assertIsInstance(midpoints2, np.ndarray)
        #self.assertIsInstance(means2, np.ndarray)
