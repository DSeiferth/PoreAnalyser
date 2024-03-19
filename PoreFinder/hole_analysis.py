import numpy as np
import MDAnalysis
import matplotlib.pyplot as plt
import pandas as pd
from MDAnalysis.analysis import hole2

#bla = '/biggin/b198/orie4254/Documents/PoreFinding_pdb/' # when loading module
import os
#bla = os.getcwd() + '/' # # when running app
#bla = os.getcwd() + '/' +'PoreFinding/' # when running locally
bla = os.path.realpath(__file__)[:-16] 
print('path hole_analysis',bla)
#bla = './' # for using streamlit app
hole_exe = bla+'hole2/hole'
print('hole_exe', hole_exe)
sph_proc = bla+'hole2/sph_process'
f_size = 18

import warnings; warnings.simplefilter('ignore')

def align_to_z(p, pdb_name, align_bool=True, sel="protein"):
    """
    rotate the principal axes of the molecule to align with Cartesian coordinate system
    """
    conf =  p + pdb_name + '.pdb'
    top = conf
    u = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')
    protein = u.select_atoms(sel)
    if align_bool:
        CA = u.select_atoms("protein and name CA")
        I = CA.moment_of_inertia()
        print('moment_of_inertia', I)
        # https://pythoninchemistry.org/ch40208/comp_chem_methods/moments_of_inertia.html
        eigenval , eigenvec = np.linalg.eig(I)
        tranform = np.linalg.inv(eigenvec)
        # rotate such that principal axis are aligned to Cartesian axis
        protein.rotate(tranform)
        # rotate 90 deg around y-axis
        R_matrix = np.zeros((3,3))
        R_matrix[0][2] = 1 
        R_matrix[1][1] = 1 
        R_matrix[2][0] = -1
        protein.rotate(R_matrix)
    protein.write(p + pdb_name + '_aligned_z.pdb')

def hole_analysis(name, path, end_radius=20, sel='protein'):
    """
    Perform hole analysis on a molecular structure and create a VMD surface for visualization.

    Parameters:
        name (str): The name of the input file.
        path (str): The path to the input file.
        end_radius (float, optional): The radius (in Angstroms) of the maximum hole to detect.
                                      Default is 20 Angstroms.
        sel (str, optional): The selection string to select the atoms to be analyzed.
                             Default is 'protein'.

    Returns:
        midpoints2 (numpy.ndarray): An array of the midpoints of the histogram bins.
        means2 (numpy.ndarray): An array of the mean values of the histogram bins.
    """
    tmpdir = path #+ 'tmpdir/'
    conf = path + name
    top = conf
    print('hole_analysis', conf)
    sys = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb') 
    ha2 = hole2.HoleAnalysis(
                                sys, select=sel,
                                cpoint='center_of_geometry',
                                executable=hole_exe,
                                tmpdir=tmpdir,
				                sph_process=sph_proc,
				                #sphpdb_file=path+name+'.sph',
                                #sphpdb=path+name+'.sph',
                                end_radius=end_radius,
				                #keep_files=True
                                cvect=[0,0,1],
				)
    ha2.run(random_seed=31415)
    try:
        #http://minium.com.au/UserGuide/stable/examples/analysis/polymers_and_membranes/hole.html
        ha2.create_vmd_surface(filename=path+name[:-4]+'.vmd') #source hole1.vmd
        print('vmd surface created for', name)
    except:
        print('ERROR: No vmd surface created for', name)

    radii2, edges2 = ha2.bin_radii(bins=100, range=None)
    means2, edges2 = ha2.histogram_radii(bins=100, range=None,
                                          aggregator=np.mean)
    midpoints2 = 0.5*(edges2[1:]+edges2[:-1])
    return midpoints2, means2

#def visualisation(name , path='/biggin/b198/orie4254/Documents/CHAP/', out=0,
#                  num_circle=24, TMD_higher=0,  TMD_lower=0, end_radius=15):
    
def analysis(names,labels, path='/biggin/b198/orie4254/Documents/CHAP/', end_radius=15,
            title='', sel='protein', legend_outside=False, 
            plot_lines=True, f_size=18, align_bool=True):
    """
    Perform hole analysis on one or more PDB files and plot the results.

    Parameters
    ----------
    names : list of str
        Names of the PDB files to analyze. If multiple files are provided, the function will
        align them to the first file in the list before analysis.
    labels : list of str
        Labels to use for the legend in the plot, corresponding to each PDB file.
    path : str, optional
        Path to the directory containing the PDB files. Default is '/biggin/b198/orie4254/Documents/CHAP/'.
    end_radius : float, optional
        End radius of the HOLE cylinder, in Angstroms. Default is 15.
    
    save : str, optional
        Name of the file to save the plot as (without extension). If not provided, the plot will not be saved.
    title : str, optional
        Title of the plot. Default is an empty string.
    sel : str, optional
        Selection string to use for the hole analysis. Default is 'protein'.
    legend_outside : bool, optional
        If True, place the legend outside of the plot. Default is False.

    align_bool: bool, optional
        If True, place align the largest principal component to z-axis. Default is True.

    Returns
    -------
    fig : matplotlib figure
        The generated plot figure.
    df : pandas DataFrame
        A DataFrame containing the results of the hole analysis, with the following columns:
        - 'Label z [A]': the z-coordinate of each point along the pore axis.
        - 'Label Radius [A]': the radius of the pore at each point.
        'Label' corresponds to the labels provided in the `labels` parameter.
    """
    
    for i in range(len(names)):
         align_to_z(p=path, pdb_name=names[i][:-4], align_bool=align_bool, sel=sel)
         names[i] = names[i][:-4] + '_aligned_z.pdb' 
    ### hole analysis ###
    aligned_path = path
    fig, ax = plt.subplots()
    plt.title(title, fontsize=f_size)
    colors = ['black', 'blue', 'orange', 'purple','green','red','gray', 'brown',
              'cyan', 'violet', 'olive', 'peru', 'slategray',
    ]
    for i in range(len(names)):
        midpoints, means = hole_analysis(name=names[i], path=aligned_path, #typ='pdb',
                                                end_radius=end_radius, sel=sel)
        rmin = min(means)
        ax.plot(midpoints, means, color=colors[i],
                label=labels[i] + r' R$_{min}$ = '+str(round(rmin,2)) + r'$\AA$')
        if i==0:
                    d = {labels[i]+' z [A]': midpoints,
                         labels[i]+' Radius [A]': means,
                        }
                    df = pd.DataFrame(data=d)
        else:
            df[labels[i]+' z [A]'] = midpoints
            df[labels[i]+' Radius [A]'] = means
    plt.ylabel(r"HOLE radius $R$ ($\AA$)", fontsize=f_size)
    plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)", fontsize=f_size)

    y1 = 0.0
    y2 = 20
    #ax.set_ylim([y1,y2])
    xlim = ax.get_xlim()
    if plot_lines:
        ax.plot(xlim, [1.15,1.15], '--',color='red') # label=r'r < 1.15 $\AA$'
        ax.plot(xlim, [2.3,2.3], '--',color='green') # label=r'r < 2.30 $\AA$'

    ax.tick_params(axis='both', which='major', labelsize=f_size)
    if legend_outside:
        plt.legend(prop={'size': f_size},  loc='upper center', bbox_to_anchor=(1.4,1.0), frameon=False)
    else:
        ax.legend(prop={'size': f_size}, loc='upper left') # loc='upper center'
        fig.tight_layout()
    #fig.savefig(path + save +'_HOLE_pathwayprofile.png', bbox_inches='tight')
    #df.to_csv(outpath + '/Pathway_HOLE_comparison_pdb.csv', sep=',', index=False, header=True)
    plt.show()
    ### visualise pathway ###
    pathways = []
    for count, name in enumerate(names):
        try:
            print(name, '### visualise pathway ###')
            print('pdbfile=path+name', path+name)
            print("sphpdb_file=path+name[:-4]+'.sph'", path+name[:-4]+".sph",)
            print(aligned_path+name[:-4]+".sph")
            ha2 =  hole2.hole(
                    pdbfile=path+name ,#+'.pdb',
                    #cpoint='center_of_geometry',
                    executable=hole_exe,
                    tmpdir=path,
                    #sph_process=sph_proc, #hole() got an unexpected keyword argument 'sph_process'
                    sphpdb_file=aligned_path+name[:-4]+".sph",
                    end_radius=end_radius,
                    keep_files=True,
                    cvect=[0,0,1],
            )
        except:
            print('ERROR with', name, 'no SPH file generated')
    return fig, df
