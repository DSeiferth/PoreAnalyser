import unittest
#import PoreFinding_pdb as pf 
import numpy as np
import sys
sys.path.append('ProbeParticleEllipsoid/')
import ellipse_lib as e_lib
#from ellipse_lib import atom, ellipse, distance_ellipse, dist_ellipse_vdwSphere, assign_radius

class atomTest(unittest.TestCase):
    """
    Test the atom class.
    """
    def test_create(self):
        """
        Tests atom creation.
        """
        a = e_lib.atom(x=1, y=2, r=4)
        self.assertIsInstance(a, e_lib.atom)
        self.assertEqual(a.r, 4)
        self.assertEqual(a.z, 0)

class ellipseTest(unittest.TestCase):
    """
    Test the ellipse class.
    """
    def test_create(self):
        e = e_lib.ellipse(a=2, b=3, cx=1, cy=2)
        self.assertIsInstance(e, e_lib.ellipse)
        self.assertEqual(e.a, 2)
        self.assertEqual(e.cz, 0)
        self.assertEqual(e.theta, 0)

    def test_on_ellipse(self):
        """
        Tests on_ellipse function.
        """
        e = e_lib.ellipse(a=2, b=3, cx=1, cy=2)
        self.assertTrue(e.on_ellipse(1, 2))
        self.assertTrue(e.on_ellipse(1, 5))
        self.assertTrue(e.on_ellipse(3, 2))
        self.assertFalse(e.on_ellipse(3, 3))

    def test_draw(self):
        """
        Tests draw function.
        """
        e = e_lib.ellipse(a=2, b=3, cx=1, cy=2)
        x, y = e.draw(res=0.01)
        self.assertEqual(len(x), len(y))
        self.assertIsInstance(x, np.ndarray)
        self.assertIsInstance(y, np.ndarray)

class distance_ellipseTest(unittest.TestCase):
    """
    Test the distance_ellipse function.
    """
    def test_distance_ellipse(self):
        """
        Tests distance_ellipse function.
        """
        semi_major = 2
        semi_minor = 3
        p = [1, 2]
        d = e_lib.distance_ellipse(semi_major, semi_minor, p)
        self.assertEqual(d[0], 1.3262624361743103)
        self.assertEqual(d[1], 2.2455094941647906)

        d = e_lib.distance_ellipse(3.0, 2.0, (1.0, 1.0))
        self.assertEqual(d[0], 1.2487110341841325)
        self.assertEqual(d[1], 1.8185123044348084)

class dist_ellipse_vdwSphereTest(unittest.TestCase):
    """
    Test the dist_ellipse_vdwSphere function.
    """
    def test_dist_ellipse_vdwSphere(self):
        """
        Tests dist_ellipse_vdwSphere function.
        """
        ellipse = e_lib.ellipse(a=3.0, b=2.0, cx=0.0, cy=0.0, theta=0.0)
        sphere = e_lib.atom(x=1.0, y=1.0, r=0.5)
        d = e_lib.dist_ellipse_vdwSphere(ellipse, sphere, plot=0)
        self.assertEqual(d, -0.5)

        ellipse = e_lib.ellipse(a=3.0, b=2.0, cx=5.0, cy=5.0, theta=0.0)
        sphere = e_lib.atom(x=1.0, y=1.0, r=0.5)
        d = e_lib.dist_ellipse_vdwSphere(ellipse, sphere, plot=0)
        self.assertEqual(d,  2.6950072040653335 ) 

class assign_radiusTest(unittest.TestCase):
    """
    Test the assign_radius function.
    """
    def test_assign_radius(self):
        """
        Tests assign_radius function.
        """
        C = e_lib.assign_radius('C')
        O = e_lib.assign_radius('O')
        S = e_lib.assign_radius('S')
        N = e_lib.assign_radius('N')
        H = e_lib.assign_radius('H')
        P = e_lib.assign_radius('P')
        self.assertEqual(C, 1.85)
        self.assertEqual(O, 1.65)
        self.assertEqual(S, 2.00)
        self.assertEqual(N, 1.75)
        self.assertEqual(H, 1.00)
        self.assertEqual(P, 2.10)