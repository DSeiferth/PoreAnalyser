import os
#bla = os.getcwd() + '/PoreFinding/'
bla = os.path.realpath(__file__)[:-15] 
print(bla)
import sys
sys.path.append(bla)
import hole_analysis as hole_analysis
from visualization import write_pdb_with_pore_surface, plt_ellipsoid_pathway, pathway_visu, st_write_ellipsoid, write_pdb_with_ellipsoid_surface, example_xy_plane, compare_volume #, render_visu
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

import conductance as conduct 

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
    Parameters for condutance estimation:
    - D_cation (float, optional): Diffusion coefficient of the cation. Default is 1.8e-9 m^2/s (value for potassium).
    - D_anion (float, optional): Diffusion coefficient of the anion. Default is 2.032e-9 m^2/s (value for chloride).
    - popt (list, optional): Parameters for the conductivity model. Default is [1.40674664, 1.25040698].
        - popt[0] (float): Scaling parameter of radius for the conductivity model (dimension 1/Angstrom).
        - popt[1] (float): Shift parameter of the sigmoid function for the conductivity model (dimenionless).
    - temp (int, optional): Temperature in Kelvin. Default is 300.
    - c_m (float, optional): Concentration in mol/l. Default is 0.15.
    - trajectory (bool, optional): Flag indicating whether the input is a trajectory. Default is False.
        - If trajectory is True, the input pdb_array should contain the path to the topology file and the trajectory file.
        - pdb_array[0]: Path to the topology file.
        - pdb_array[1]: Path to the trajectory file.
    - traj_frames (int, optional): Number of frames to extract from the trajectory. Default is 1.

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
    - popt (list): Parameters for the conductivity model.
    - bulk conductivity (float): Bulk conductivity of the system.

    Methods:
    - hole_analysis: Perform hole analysis on the set of PDB models.
    - ellipsoid_analysis: Perform ellipsoid analysis on a specific PDB model.
    - plt_pathway_ellipsoid: Plot ellipsoid analysis results for a specific model.
    - pathway_visualisation: Visualize the pathway for a specific model.
    - conductance_estimation: Estimate the conductance of the pore using a conductivity model.
    - plt_trajectory_average: Plot the trajectory average of the radius / radii profile.

    Example:
    >>> pdb_models = ['model1.pdb', 'model2.pdb']
    >>> pore_analysis = PoreAnalysis(pdb_array=pdb_models, opt_method='nelder-mead')
    >>> pore_analysis.hole_analysis()
    >>> pore_analysis.ellipsoid_analysis(index_model=0)
    >>> pore_analysis.plt_pathway_ellipsoid(index_model=0)
    >>> pore_analysis.pathway_visualisation(index_model=0)
    >>> pore_analysis.pathway_rendering(index_model=0)
    >>> pore_analysis.conductance_estimation(index_model=0)

    Example with trajectory:
    >>> pdb_models = [fname+'.tpr', fname+'.xtc']
    >>> pore_analysis = PoreAnalysis(pdb_array=pdb_models, trajectory=True, traj_frames=10)
    >>> pore_analysis.hole_analysis()
    >>> pore_analysis.plt_trajectory_average(HOLE_profile=True)
    >>> for i in range(10): pore_analysis.ellipsoid_analysis(index_model=i)
    >>> pore_analysis.plt_trajectory_average(HOLE_profile=False)
    """
    def __init__(self, pdb_array, opt_method='nelder-mead',
                 align_bool=True, end_radius=15, pathway_sel='protein',
                 path_save = '', 
                 num_circle=24, clipping=100,
                 D_cation=1.8e-9, D_anion=2.032e-9,
                 popt = [1.40674664, 1.25040698], 
                 temp=300, c_m=0.15,
                 trajectory=False,
                 traj_frames=1, 
                 ):
        self.trajectory = trajectory
        if trajectory:
            self.traj_frames = traj_frames 
            # write out frames from trajectory file
            top = pdb_array[0]
            conf = pdb_array[1]
            fname = conf.split('/')[-1].split('.')[0]
            u = MDAnalysis.Universe(top, conf)
            protein = u.select_atoms(pathway_sel)
            indices = np.arange(0, len(u.trajectory), int((len(u.trajectory))/(traj_frames)))[:traj_frames]
            self.pdb_array = []
            for i in indices:
                u.trajectory[i]
                s = path_save+fname+str(i)+'.pdb'
                protein.write(s)
                self.pdb_array.append(s)
        else:
            self.pdb_array = pdb_array
        print('pdb_array, self'  , self.pdb_array, 'input',  pdb_array )    
        self.align_bool = align_bool
        self.end_radius = end_radius
        self.pathway_sel = pathway_sel
        self.path_save = path_save
        labels = []
        names_aligned = []
        for ele in self.pdb_array:
            splits = ele.split('/')[-1]
            labels.append(splits.split('.')[0] ) #splits[:-4]
            names_aligned.append(ele.split('.')[-2] +'_aligned_z.pdb') # ele[:-4]
        self.labels = labels
        self.names_aligned = names_aligned
        print('self.labels', self.labels )
        print('self.names_aligned', self.names_aligned)

        self.opt_method = opt_method

        self.num_circle = num_circle  # for point cloud visualization of pore surface
        self.clipping = clipping # 3d visualisation of pore surface

        self.hole_fig = None 
        self.hole_df = None 

        self.ellipsoid_dfs = {}


        ### conductance estimation ###
        self.popt = popt
        e = 1.6022*1e-19 # C
        J2kcal = 1/4184 # 1 kcalth = 4184 J
        kT = 1.380*1e-23 *temp #* J2kcal# J/K
        # mobility
        mu_Na = D_cation/kT
        mu_Cl = D_anion/kT
        # concentration
        l2m3 = 0.001
        mol = 6.022*1e23 # particles
        c = c_m *1/l2m3 *mol # mol/ m^3
        self.bulk_conductivity = e*e*(mu_Na+mu_Cl)*c
        print('bulk_conductivity', self.bulk_conductivity)

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
        Figure  and dataframe

        Notes:
        This method performs hole analysis on the set of PDB models and generates visualizations.
        The results are stored in the attributes hole_fig and hole_df.

        Example:
        >>> pore_analysis = PoreAnalysis(pdb_array=['model1.pdb', 'model2.pdb'])
        >>> pore_analysis.hole_analysis()
        """
        fig, df = hole_analysis.analysis(self.pdb_array, labels=self.labels, 
                                          path=self.path_save, #'', 
                                          end_radius=self.end_radius, 
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
        return fig, df  
    
    def plt_trajectory_average(self, num_bins=100, f_size=20, title='', HOLE_profile=True):
        """
        Plot the trajectory average of the hole radius.
        Parameters:
        - num_bins (int, optional): Number of bins for the plot. Default is 100.
        - f_size (int, optional): Font size for the plot. Default is 20.
        - title (str, optional): Title for the plot. Default is an empty string.
        - HOLE_profile (bool, optional): Flag indicating whether to plot the HOLE 
            profile or the PoreAnalysor profile. Default is True.
        Returns:
        Figure  and dataframe
        """
        z = np.array([])
        if HOLE_profile:
            r = np.array([])
            for l in self.labels:
                z = np.append(z, self.hole_df[l+' z [A]'])
                r = np.append(r, self.hole_df[l+' Radius [A]'])
            df = pd.DataFrame({'z':z, 'r':r})
            bin_edges = pd.cut(df['z'], bins=num_bins)
            grouped = df.groupby(bin_edges)
            result = grouped.agg({
                'z': ['mean', 'std'],
                'r': ['mean', 'std'], })
            self.av_hole = result 
        else:
            a = np.array([])
            b = np.array([])
            for key in self.ellipsoid_dfs:
                df_PA = self.ellipsoid_dfs[key]
                z = np.append(z, df_PA.z)
                a = np.append(a, df_PA.a)
                b = np.append(b, df_PA.b)
            df = pd.DataFrame({'z':z, 'a':a, 'b':b})
            bin_edges = pd.cut(df['z'], bins=num_bins)
            grouped = df.groupby(bin_edges)
            result = grouped.agg({
                'z': ['mean', 'std'],
                'a': ['mean', 'std'],
                'b': ['mean', 'std'],
            })
            self.av_PA = result 
            
        fig, ax = plt.subplots()
        ax.set_title(title+' trajectory average', fontsize=f_size)
        if HOLE_profile:
            plt.plot(result['z']['mean'], result['r']['mean'], label='r', color='blue'  )
            plt.fill_between(
                result['z']['mean'], result['r']['mean']-result['r']['std'],  result['r']['mean']+result['r']['std'],
                color='blue', alpha=0.2,
            )
        else:
            plt.plot(result['z']['mean'], result['a']['mean'], label='a', color='blue'  )
            plt.plot(result['z']['mean'], result['b']['mean'], label='b', color='orange'  )
            plt.fill_between(
                result['z']['mean'], result['a']['mean']-result['a']['std'],  result['a']['mean']+result['a']['std'],
                color='blue', alpha=0.2,
            )
            plt.fill_between(
                result['z']['mean'], result['b']['mean']-result['b']['std'],  result['b']['mean']+result['b']['std'],
                color='orange', alpha=0.2,
            )
        ax.set_ylim([0, self.end_radius+10])
        ax.set_ylabel('Radius ($\AA$)', fontsize=f_size)
        ax.set_xlabel('z ($\AA$)', fontsize=f_size)
        ax.tick_params(axis='both', which='major', labelsize=f_size)
        ax.legend(prop={'size': f_size})
        fig.tight_layout()
        plt.show()
        return fig, result
    
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
        datframe df_res

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
                                     fname=self.names_aligned[index_model]+'_pathway_ellipse.txt', num_circle = self.num_circle)
        return df_res

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
        xyzview = pathway_visu(path=self.path_save, #path='', 
                               name=self.names_aligned[index_model], 
                               f_end=f_end, pathway_sel=self.pathway_sel,
                               clipping=self.clipping,
                               )
        # f_end='_ellipsoid.pdb' f_end='_circle.pdb'
        return xyzview
        
    def conductance_estimation(self, index_model=0, f_size=15):
        """
        Estimate the conductance of the pore using a conductivity model.
        Parameters:
        - index_model (int, optional): Index of the model in the pdb_array. Default is 0.
        - f_size (int, optional): Font size for the plot. Default is 15.

        Returns:
        tuple: Tuple containing the conductance values in pS for the pore.
        1. conductance with bulk conductivity and spherical probe particle (hole)
        2. conductance with bulk conductivity and ellipsoid probe particle (PoreAnalyser)
        3. conductance with conductivity model and spherical probe particle (hole)
        4. conductance with conductivity model and ellipsoid probe particle (PoreAnalyser)
        5. fig : matplotlib.figure.Figure: Figure object for the plot (resistance vs z)
        """
        pS = 1e-12
        system = self.labels[index_model]
        print(system)
        res1 = np.loadtxt(self.path_save+system+'_aligned_z.pdb_pathway_ellipse.txt', 
                                comments='#', delimiter=',')
        df_res1 = pd.DataFrame(data=res1, columns=['x', 'y', 'z', 'a', 'b', 'theta'])
        df_res1.sort_values('z', inplace=True)

        hole1, R = conduct.bullk_conduct(z=df_res1['z'],a=df_res1['b'], b=df_res1['b'],conduct=self.bulk_conductivity)
        pf1, R_pf_bulk = conduct.bullk_conduct(z=df_res1['z'],a=df_res1['a'],b=df_res1['b'],conduct=self.bulk_conductivity)
        hole_c, R_hole_c, facs_hole = conduct.no_bulk_conduct(z=df_res1['z'],a=df_res1['b'],b=df_res1['b'], plot=False, popt=self.popt, conduct=self.bulk_conductivity) 
        pf1_c, R_pf_c, facs_pf = conduct.no_bulk_conduct(z=df_res1.z,a=df_res1.a,b=df_res1.b, plot=False, popt=self.popt, conduct=self.bulk_conductivity)

        print('conductance (pS)','hole', hole1/pS, 'pf', pf1/pS,)
        print('conductance (pS)', 'hole_c', hole_c/pS, 'pf_c', pf1_c/pS,)
        print() 
        fig, ax = plt.subplots(nrows=1, ncols=2, 
                       sharex=True,
                       #sharey='row', 
                       figsize=(15, 10))
        ax[0].plot(df_res1['z'][1:], facs_hole, label=r'HOLE with conductivity model')
        ax[0].plot(df_res1['z'][1:], facs_pf, label=r'PoreAnalyser with conductivity model')

        Ohm_2GOhm = 1e-9
        ax[1].plot(df_res1['z'][1:], np.array(R_hole_c)*Ohm_2GOhm, label=r'HOLE with conductivity model')
        ax[1].plot(df_res1['z'][1:], np.array(R_pf_c)*Ohm_2GOhm, label=r'PoreAnalyser with conductivity model')
        ax[1].plot(df_res1['z'][1:], np.array(R)*Ohm_2GOhm, label=r'HOLE with $\kappa_{bulk}$')
        ax[1].plot(df_res1['z'][1:], np.array(R_pf_bulk)*Ohm_2GOhm, label=r'PoreAnalyser with $\kappa_{bulk}$')
        ax[1].legend(prop={'size': f_size*0.9}) 

        ax[0].set_title('Conductivity along pore axis', fontsize=f_size, loc='center')
        ax[1].set_title('Resistance along pore axis', fontsize=f_size, loc='center')
        ax[0].set_xlabel("z ($\AA$)", fontsize=f_size)
        ax[0].set_ylabel(r"$\kappa$(a,b)/\kappa_{bulk}$", fontsize=f_size)
        ax[1].set_xlabel("z ($\AA$)", fontsize=f_size)
        ax[1].set_ylabel(r"Resistance ($G\Omega$)", fontsize=f_size)
        ax[0].tick_params(axis='both', which='major', labelsize=f_size)
        ax[1].tick_params(axis='both', which='major', labelsize=f_size)
        #fig.tight_layout()
        plt.show()
        return hole1/pS, pf1/pS, hole_c/pS, pf1_c/pS, fig 
