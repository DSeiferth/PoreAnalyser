import unittest
import numpy as np
import sys
import MDAnalysis
import pandas as pd
sys.path.append('ProbeParticleEllipsoid/')
import ellipsoid_optimisation as e_opt 
import ellipse_lib as e_lib
import os
bla = os.getcwd() + '/' +'PoreAnalyser/' # when running locally


class penalty_overlap_4dim_Test(unittest.TestCase):
    """
    Test the penalty_overlap_4dim function.
    """
    def test_penalty_overlap_4dim(self):
        """
        Tests penalty_overlap_4dim function.
        """
        # Create test data
        a_vec = [e_lib.atom(x=3.8, y=4.8, r=2), e_lib.atom(-1, -3, r=2.7), 
                 e_lib.atom(-4, 4, r=1), e_lib.atom(7, 2, r=1), 
                 e_lib.atom(7.9, -1, r=1.5), e_lib.atom(-3, 1, r=1), 
                 e_lib.atom(-0.0, 4.8, r=1.3), e_lib.atom(3, -2.9, r=1.0)]
        probe = e_lib.atom(2, 1, r=2)
        # Test function for no overlap
        x = [probe.r, 0, probe.x, probe.y]
        args = [probe.r, a_vec, True]
        penalty = e_opt.penalty_overlap_4dim(x, args)
        # Test results for no overlap
        self.assertEqual(penalty, -probe.r)  # no  overlap => penalty is negative radius

        # ## overlap with stop_Loop ###
        x = [3, 0, probe.x, probe.y]
        args = [probe.r, a_vec, True]
        penalty = e_opt.penalty_overlap_4dim(x, args)
        self.assertEqual(penalty, 1000000000.0)

        # ## overlap with stop_Loop=False  ###
        x = [3, 0, probe.x, probe.y]
        args = [probe.r, a_vec, False]
        penalty = e_opt.penalty_overlap_4dim(x, args)
        self.assertEqual(penalty, 0.01671842700018855)


class neighbor_vec_Test(unittest.TestCase):
    def test_neighbor_vec(self):
        """
        Tests neighbor_vec function.
        """
        # ## load data from GlyR model ###
        p = bla+'pdb_models/'
        name = 'test_7tu9_aligned_z'
        end_radius = 20
        
        conf = p + name + '.sph'
        top = conf
        sph = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

        conf = p + name + '.pdb' 
        top = conf
        u = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

        mer = MDAnalysis.Merge(u.atoms, sph.atoms).select_atoms('name *')

        radii = sph.atoms.occupancies
        resids = sph.atoms.resids
        x_coordinates = []
        y_coordinates = []
        z_coordinates = []
        for i, pos in enumerate(sph.atoms.positions):
            x_coordinates.append(pos[0])
            y_coordinates.append(pos[1])
            z_coordinates.append(pos[2])

        x_coordinates = np.array(x_coordinates)
        y_coordinates = np.array(y_coordinates)
        z_coordinates = np.array(z_coordinates)

        ind_sort = np.argsort(z_coordinates)
        x_coordinates = x_coordinates[ind_sort]
        y_coordinates = y_coordinates[ind_sort]
        z_coordinates = z_coordinates[ind_sort]
        radii = radii[ind_sort]
        resids = resids[ind_sort]

        d = {'x': x_coordinates, 'y': y_coordinates, 'z': z_coordinates, 'r': radii, 'resid': resids}
        df = pd.DataFrame(data=d)
        df2 = df[(df['r'] < end_radius)]

        # ## prepare input for function ###
        dataframe = df2
        index = 500
        universe = mer

        rID = dataframe['resid'].loc[index]
        probe1 = universe.select_atoms('resid '+str(rID) + ' and resname SPH')

        probe = e_lib.atom(dataframe['x'].loc[index],dataframe['y'].loc[index],
                        z=dataframe['z'].loc[index],
                        r=dataframe['r'].loc[index])
            
        n_xy_fac = 1.6
        out = 0 
        pathway_sel = 'protein'

        a_vec, neighbour_labels, n_xy = e_opt.neighbor_vec(universe, probe, probe1, n_xy_fac, 
                                                        out=out, pathway_sel=pathway_sel)
        
        self.assertEqual(len(a_vec), 50)
        self.assertEqual(len(neighbour_labels), 50)
        self.assertEqual(neighbour_labels[0], 'SER 294')
        self.assertEqual(n_xy, 10.15625 )

        self.assertEqual(a_vec[0].r, 1.1284209732286377 )
        self.assertEqual(a_vec[10].r, 1.310291834203602 )
        self.assertEqual(a_vec[40].r, 0.3332456537007631  )
