import unittest
import numpy as np
import pandas as pd 
import conductance as conduct 
import os
# bla = '' # when running locally 
bla = os.getcwd() + '/'  +'PoreAnalyser/'

popt = [1.40674664, 1.25040698]
pS = 1e-12

class conductanceTest(unittest.TestCase):
    """
    Test the conductance functions.
    """
    def test_bullk_conduct(self):
        """
        Tests bullk_conduct function 
        """
        hole1_vec = [77.83958791322752, 151.6641306537728, 127.7991033324672 ]
        pf1_vec = [100.90087370272876, 227.90060005260435, 172.62682152696942 ]
        for count, system in enumerate(['7tu9', '7tvi', '8fe1']):
            res1 = np.loadtxt(bla+'pdb_models/'+system+'_aligned_z.pdb_pathway_ellipse.txt', 
                                comments='#', delimiter=',')
            df_res1 = pd.DataFrame(data=res1, columns=['x', 'y', 'z', 'a', 'b', 'theta'])
            df_res1.sort_values('z', inplace=True)

            hole1, R = conduct.bullk_conduct(z=df_res1['z'],a=df_res1['b'], b=df_res1['b'])
            pf1, R_pf_bulk = conduct.bullk_conduct(z=df_res1['z'],a=df_res1['a'],b=df_res1['b'])
        
            self.assertEqual(hole1/pS, hole1_vec[count])
            self.assertEqual(pf1/pS, pf1_vec[count])

    def test_no_bulk_conduct(self):
        """
        Tests no_bulk_conduct function 
        """
        hole_vec = [9.921242336202274, 23.180663879581154, 17.79738713224399 ]
        pf_vec = [13.478509637919807, 41.066684695303124, 25.708421014160677 ]
        for count, system in enumerate( ['7tu9', '7tvi', '8fe1']):
            res1 = np.loadtxt(bla+'pdb_models/'+system+'_aligned_z.pdb_pathway_ellipse.txt', 
                                comments='#', delimiter=',')
            df_res1 = pd.DataFrame(data=res1, columns=['x', 'y', 'z', 'a', 'b', 'theta'])
            df_res1.sort_values('z', inplace=True)
            
            hole_c, R_hole_c, facs_hole = conduct.no_bulk_conduct(z=df_res1['z'],a=df_res1['b'],b=df_res1['b'], plot=False, popt=popt) 
            pf1_c, R_pf_c, facs_pf = conduct.no_bulk_conduct(z=df_res1.z,a=df_res1.a,b=df_res1.b, plot=False, popt=popt)
        
            self.assertEqual(hole_c/pS, hole_vec[count])
            self.assertEqual(pf1_c/pS, pf_vec[count])
