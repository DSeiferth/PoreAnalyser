import os
#bla = os.getcwd() + '/PoreFinding/'
bla = os.path.realpath(__file__)[:-14] 
print(bla)
import sys
sys.path.append(bla)
import hole_analysis as hole_analysis
from visualization import write_pdb_with_pore_surface, plt_ellipsoid_pathway, pathway_visu, st_write_ellipsoid, write_pdb_with_ellipsoid_surface, example_xy_plane, compare_volume, render_visu
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
    """
    Class for pore analysis of a set of pdb models.
    Parameters:
    - pdb_array (list): List of file paths to the input PDB models.
    - opt_method (str, optional): Optimization method for ellipsoid fitting. Default is 'nelder-mead'.
    - align_bool (bool, optional): Flag indicating whether to align the models. Default is True.
    - end_radius (float, optional): Radius of the spherical probe particle. Default is 15.
    - pathway_sel (str, optional): Atom selection string for identifying the pathway. Default is 'protein'.
    - path_save (str, optional): Path to save analysis results. Default is ''.
    - num_circle (int, optional): Number of circles for point cloud visualization of pore surface. Default is 24.
    - clipping (int, optional): Clipping value for 3D visualization of pore surface. Default is 100.

    Attributes:
    - pdb_array (list): List of file paths to the input PDB models.
    - align_bool (bool): Flag indicating whether to align the models.
    - end_radius (float): Radius of the spherical probe particle.
    - pathway_sel (str): Atom selection string for identifying the pathway.
    - path_save (str): Path to save analysis results.
    - labels (list): List of labels derived from the input PDB file paths.
    - names_aligned (list): List of aligned PDB file paths.
    - opt_method (str): Optimization method for ellipsoid fitting.
    - num_circle (int): Number of circles for point cloud visualization of pore surface.
    - clipping (int): Clipping value for 3D visualization of pore surface.
    - hole_fig (matplotlib.figure.Figure): Figure object for the hole analysis.
    - hole_df (pandas.DataFrame): DataFrame containing hole analysis results.
    - ellipsoid_dfs (dict): Dictionary containing ellipsoid analysis results for each model.

    Methods:
    - hole_analysis: Perform hole analysis on the set of PDB models.
    - ellipsoid_analysis: Perform ellipsoid analysis on a specific PDB model.
    - plt_pathway_ellipsoid: Plot ellipsoid analysis results for a specific model.
    - pathway_visualisation: Visualize the pathway for a specific model.
    - pathway_rendering: Render the pathway for a specific model.

    Example:
    >>> pdb_models = ['model1.pdb', 'model2.pdb']
    >>> pore_analysis = PoreAnalysis(pdb_array=pdb_models, opt_method='nelder-mead')
    >>> pore_analysis.hole_analysis()
    >>> pore_analysis.ellipsoid_analysis(index_model=0)
    >>> pore_analysis.plt_pathway_ellipsoid(index_model=0)
    >>> pore_analysis.pathway_visualisation(index_model=0)
    >>> pore_analysis.pathway_rendering(index_model=0)
    """
    def __init__(self, pdb_array, opt_method='nelder-mead',
                 align_bool=True, end_radius=15, pathway_sel='protein',
                 path_save = '', 
                 num_circle=24, clipping=100,
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

        self.num_circle = num_circle  # for point cloud visualization of pore surface
        self.clipping = clipping # 3d visualisation of pore surface

        self.hole_fig = None 
        self.hole_df = None 

        self.ellipsoid_dfs = {}

    def hole_analysis(self, plot_lines=True, legend_outside=False, title='', f_size=15):
        """
        Perform hole analysis on the set of PDB models.
        HOLE uses a spherical probe particle.

        Parameters:
        - plot_lines (bool, optional): Flag indicating whether to plot lines. Default is True.
        - legend_outside (bool, optional): Flag indicating whether to place the legend outside the plot. Default is False.
        - title (str, optional): Title for the plot. Default is an empty string.
        - f_size (int, optional): Font size for the plot. Default is 15.

        Returns:
        None

        Notes:
        This method performs hole analysis on the set of PDB models and generates visualizations.
        The results are stored in the attributes hole_fig and hole_df.

        Example:
        >>> pore_analysis = PoreAnalysis(pdb_array=['model1.pdb', 'model2.pdb'])
        >>> pore_analysis.hole_analysis()
        """
        fig, df = hole_analysis.analysis(self.pdb_array, labels=self.labels, 
                                          path='', end_radius=self.end_radius, 
                                          title=title,
                                          legend_outside=legend_outside, plot_lines=plot_lines, 
                                          f_size=f_size, align_bool=self.align_bool, 
                                          sel=self.pathway_sel
                                            )
        self.hole_fig = fig
        self.hole_df = df 
        for i in range(len(self.pdb_array)):
            write_pdb_with_pore_surface(path=self.path_save, name=self.names_aligned[i], 
                                        end_radius=self.end_radius, num_circle = self.num_circle)
    
    def ellipsoid_analysis(self, index_model=0, 
                           plot_lines=True, legend_outside=False, title='', f_size=15):
        """
        Perform ellipsoid analysis on a specific PDB model.

        Parameters:
        - index_model (int, optional): Index of the model in the pdb_array. Default is 0.
        - plot_lines (bool, optional): Flag indicating whether to plot lines. Default is True.
        - legend_outside (bool, optional): Flag indicating whether to place the legend outside the plot. Default is False.
        - title (str, optional): Title for the plot. Default is an empty string.
        - f_size (int, optional): Font size for the plot. Default is 15.

        Returns:
        None

        Notes:
        This method performs ellipsoid analysis on a specific PDB model and generates visualizations.
        The results are stored in the ellipsoid_dfs attribute.

        Example:
        >>> pore_analysis = PoreAnalysis(pdb_array=['model1.pdb', 'model2.pdb'])
        >>> pore_analysis.ellipsoid_analysis(index_model=0)
        """
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

        self.ellipsoid_dfs[self.labels[index_model]] = df_res

        write_pdb_with_ellipsoid_surface(p='', pdbname=self.names_aligned[index_model], 
                                     fname=self.names_aligned[0]+'_pathway_ellipse.txt', num_circle = self.num_circle)

    def plt_pathway_ellipsoid(self, index_model=0, title='', f_size=15):
        """
        Plot ellipsoid analysis results for a specific model.

        Parameters:
        - index_model (int, optional): Index of the model in the pdb_array. Default is 0.
        - title (str, optional): Title for the plot. Default is an empty string.
        - f_size (int, optional): Font size for the plot. Default is 15.

        Returns:
        matplotlib.figure.Figure: Figure object for the plot.

        Example:
        >>> pore_analysis = PoreAnalysis(pdb_array=['model1.pdb', 'model2.pdb'])
        >>> pore_analysis.ellipsoid_analysis(index_model=0)
        >>> fig = pore_analysis.plt_pathway_ellipsoid(index_model=0)
        """
        df_res = self.ellipsoid_dfs[self.labels[index_model]]
        fig = plt_ellipsoid_pathway(df_res, f_size=f_size, title=title, end_radius=self.end_radius)
        return fig 
    
    def pathway_visualisation(self, index_model=0, f_end='_circle.pdb'):
        """
        Visualize the pathway for a specific model.

        Parameters:
        - index_model (int, optional): Index of the model in the pdb_array. Default is 0.
        - f_end (str, optional): File ending for the visualization file. Default is '_circle.pdb'.

        Returns:
        xyzview: Pathway visualization object.

        Example:
        >>> pore_analysis = PoreAnalysis(pdb_array=['model1.pdb', 'model2.pdb'])
        >>> xyzview = pore_analysis.pathway_visualisation(index_model=0)
        """
        xyzview = pathway_visu(path='', name=self.names_aligned[index_model], 
                               f_end=f_end, pathway_sel=self.pathway_sel,
                               clipping=self.clipping,
                               )
        # f_end='_ellipsoid.pdb' f_end='_circle.pdb'
        return xyzview
    
    def pathway_rendering(self, index_model=0, f_end='_circle.pdb', outname='out'):
        """
        Render the pathway for a specific model.

        Parameters:
        - index_model (int, optional): Index of the model in the pdb_array. Default is 0.
        - f_end (str, optional): File ending for the visualization file. Default is '_circle.pdb'.
        - outname
        """
        render_visu(path='', name=self.names_aligned[index_model], 
                    f_end=f_end, outname=outname, streamlit=False)
