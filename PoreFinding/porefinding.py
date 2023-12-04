import os
bla = os.getcwd() + '/PoreFinding/'
print(bla)
import sys
sys.path.append(bla)
import hole_analysis as hole_analysis
import MDAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    import multiprocessing 
    print("Number of processors: ", multiprocessing.cpu_count())
    parallel = True
except:
    parallel = False
    print("Could not import multiprocessing library, => multiprocessing disabled")
sys.path.append(bla+'ProbeParticleEllipsoid/')
from ellipsoid_optimisation import ellipsoid_pathway

class PoreAnalysis():
    def __init__(self, pdb_array, opt_method='nelder-mead',
                 align_bool=True, end_radius=15, pathway_sel='protein',
                 path_save = '', 
                 ):
        self.pdb_array = pdb_array
        self.align_bool = align_bool
        self.end_radius = end_radius
        self.pathway_sel = pathway_sel
        self.path_save = path_save
        labels = []
        names_aligned = []
        for ele in pdb_array:
            splits = ele.split('/')[-1]
            labels.append(splits[:-4])
            names_aligned.append(ele[:-4]+'_aligned_z.pdb')
        self.labels = labels
        self.names_aligned = names_aligned

        self.opt_method = opt_method

        self.hole_fig = None 
        self.hole_df = None 

        self.ellipsoid_dfs = {}

    def hole_analysis(self, plot_lines=True, legend_outside=False, title='', f_size=15):
        """
        spherical probe particle
        """
        fig , df = hole_analysis.analysis(self.pdb_array, labels=self.labels, 
                                          path='', end_radius=self.end_radius, 
                                          title=title,
                                          legend_outside=legend_outside, plot_lines=plot_lines, 
                                          f_size=f_size, align_bool=self.align_bool, 
                                          sel=self.pathway_sel
                                            )
        self.hole_fig = fig
        self.hole_df = df 
    
    def ellipsoid_analysis(self, index_model=0, 
                           plot_lines=True, legend_outside=False, title='', f_size=15):
        ellipsoid_pathway(p=self.path_save, 
                        pdb_name = self.names_aligned[index_model], 
                        sph_name = self.names_aligned[index_model][:-4], 
                        slice_dz=4, parallel=parallel,  
                        num_processes=None, timeout=6, 
                        start_index = 1, end_radius=self.end_radius-1,
                        out = 0,
                        n_xy_fac = 3,#1.6,
                        pathway_sel=self.pathway_sel,
                        opt_method=self.opt_method,
                    )
        res = np.loadtxt(self.path_save + self.names_aligned[index_model]+ '_pathway_ellipse.txt', 
                    comments='#', delimiter=',')
        df_res = pd.DataFrame(data=res, columns=['x', 'y', 'z', 'a', 'b', 'theta'])
        df_res.sort_values('z', inplace=True)

        self.ellipsoid_dfs[self.names_aligned[index_model][:-4]] = df_res
