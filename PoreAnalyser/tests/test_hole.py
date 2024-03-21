import unittest
import hole_analysis as hole_analysis
import numpy as np
import os
bla = os.getcwd() + '/'  +'PoreAnalyser/' # when running locally


class HoleTest(unittest.TestCase):
    """
    Tests the hole_analysis function.
    """
    def test_hole_analysis(self):
        """
        Tests hole_analysis function.
        """
        path = bla+'pdb_models/'
        name = '7tvi.pdb'
        midpoints2, means2 = hole_analysis.hole_analysis(name, path, end_radius=20, sel='protein')

        self.assertEqual(len(midpoints2), len(means2))
        self.assertIsInstance(midpoints2, np.ndarray)
        self.assertIsInstance(means2, np.ndarray)

# OSError: [Errno 8] Exec format error: '/home/runner/work/PoreFinding_pdb/PoreFinding_pdb/hole2/hole'
# File "/home/runner/work/PoreFinding_pdb/PoreFinding_pdb/tests/test_hole.py", line 18, in test_hole_analysis
#    midpoints2, means2 = hole_analysis.hole_analysis(name, path, end_radius=20, sel='protein')

# locally: file hole2/hole
# hole2/hole: ELF 64-bit LSB executable, x86-64, version 1 (GNU/Linux), statically linked, for GNU/Linux 2.6.24, BuildID[sha1]=b539a1449b8ef241cbfe3c560c6c25963eea2725, not stripped
# on github workflow:
# 
