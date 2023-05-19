import numpy as np
import MDAnalysis
import nglview as nv
import matplotlib.pyplot as plt
import pandas as pd
from MDAnalysis.analysis import hole2


hole_exe = 'hole2//hole'
sph_proc = 'hole2/sph_process'
f_size = 18

import warnings; warnings.simplefilter('ignore')

from MDAnalysis.analysis import align
def align_and_write(universe, reference, names, out_path,
                   TMD_lower=0, TMD_higher=0):
    """
    Aligns the `universe` structure to the `reference` structure using the `alignto` method of `MDAnalysis`.
    Writes the aligned structures to `out_path` with the names specified in `names`.
    
    Args:
    universe (MDAnalysis Universe object): The structure to be aligned.
    reference (MDAnalysis Universe object): The reference structure to align to.
    names (list of str): A list of two strings, where the first string represents the name of the aligned structure and 
                         the second string represents the name of the reference structure.
                         names=[names[i] ,names[0]]: 2nd name refers to reference
    out_path (str): The path where the output files will be written.
    TMD_lower (int, optional): The lower bound of the transmembrane domain of the protein structure. Default is 0.
    TMD_higher (int, optional): The upper bound of the transmembrane domain of the protein structure. Default is 0.
    
    Returns:
    tuple: A tuple of `names` and the second element of `rmsds`.
    """
    try:
        CA_2align = universe.select_atoms('name CA ') #and chainID A and segid seg_0_PROA
        CA_ref = reference.select_atoms('name CA ') # and chainID A
        print('calphas model 1, reference (model 0)', len(CA_2align), len(CA_ref))
        align_IDs = np.unique(CA_2align.resids)
        ref_IDs = np.unique(CA_ref.resids)
        mask1 = np.isin(align_IDs, ref_IDs)
        mask2 = np.isin(ref_IDs, align_IDs)
        chains = '(chainID A )' # for glycine
        print(np.unique(align_IDs[mask1] == ref_IDs[mask2]),
              'number of common c-alphas', len(ref_IDs[mask2]),  len(align_IDs[mask1]))
        select = 'name CA and ('
        for i in ref_IDs[mask2]:
            select = select + 'resid '+str(i) + ' or '
        select = select + 'resid '+str(ref_IDs[mask2][0])+' ) ' #'and ' + chains
        sel1 = universe.select_atoms(select)
        sel2 = reference.select_atoms(select)
        print('selection in reference', len(sel2), 'selection in structure to be aligned',
              len(sel1))

        rmsds = align.alignto(universe.atoms,  # mobile
                          reference.atoms,  # reference
                          select=select, # selection to operate on
                          #match_atoms=True # whether to match atoms
                         )
        print('returns (old_rmsd, new_rmsd)',rmsds)

        sel1 = 'protein and prop  z<'+str(TMD_higher)
        sel2 = ' and prop z>'+str(TMD_lower)

        aligned = MDAnalysis.Merge(universe.atoms).select_atoms('protein')
        if TMD_higher!=0 and TMD_lower!=0 :
            TM = aligned.select_atoms(sel1+sel2)
            TM.write(out_path+  names[0]+'.pdb')
        else:
            aligned.atoms.write(out_path+  names[0]+'.pdb')

        ref = MDAnalysis.Merge(reference.atoms).select_atoms('protein')
        if TMD_higher!=0 and TMD_lower!=0 :
            TM = ref.select_atoms(sel1+sel2)
            TM.write(out_path+  names[1]+'.pdb')
        else:
            ref.atoms.write(out_path+names[1]+'.pdb')
        print('aligned', names[0], 'with', names[1], 'as reference')
    except:
        print('ERROR: Alignment did not work')
        print('usually the following error:')
        print('electionError: Reference and trajectory atom selections do not contain the same number of atoms: ')
        print('Printing the orignal models as pdb')
        universe.atoms.write(out_path+  names[0]+'.pdb')
        reference.atoms.write(out_path+names[1]+'.pdb')
        return names, [0,0,0]

    return names, rmsds[1]

def visualise_aligned(names, path):
    '''
    Args:
    names (list of str): A list of two strings, where the first string represents the name of the aligned structure and 
                         the second string represents the name of the reference structure.
    path 
    Returns:
    a nglview object
    '''
    conf =  path +  names[0] + '.pdb'
    top = conf
    u = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')
    conf =  path +  names[1] + '.pdb'
    top = conf
    ref = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

    mer = MDAnalysis.Merge(u.atoms, ref.atoms).select_atoms('protein') #name *
    v = nv.show_mdanalysis(mer.atoms)
    return v

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
    sys = MDAnalysis.Universe(top, conf) 
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
				)
    ha2.run(random_seed=31415)
    try:
        #http://minium.com.au/UserGuide/stable/examples/analysis/polymers_and_membranes/hole.html
        ha2.create_vmd_surface(filename=path+name+'.vmd') #source hole1.vmd
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
            TMD_lower=0, TMD_higher=0, save='', title='', sel='protein', legend_outside=False, 
            plot_lines=True, f_size=18):
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
    TMD_lower : float, optional
        Lower bound for the TMD selection, in Angstroms along the z-axis. If both TMD_lower and TMD_higher
        are non-zero, only atoms within the TMD range will be selected. Default is 0.
    TMD_higher : float, optional
        Upper bound for the TMD selection, in Angstroms along the z-axis. If both TMD_lower and TMD_higher
        are non-zero, only atoms within the TMD range will be selected. Default is 0.
    save : str, optional
        Name of the file to save the plot as (without extension). If not provided, the plot will not be saved.
    title : str, optional
        Title of the plot. Default is an empty string.
    sel : str, optional
        Selection string to use for the hole analysis. Default is 'protein'.
    legend_outside : bool, optional
        If True, place the legend outside of the plot. Default is False.

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
    ### align model ###
    conf =  path + names[0] 
    print('conf', conf)
    #! ls $conf
    top = conf
    ref = MDAnalysis.Universe(top, conf,) 

    if len(names) == 1 :
        print('ref', len(ref.atoms),ref)
        if TMD_lower!=0 and TMD_higher!=0:
            sel1 = 'protein and prop  z<'+str(TMD_higher)
            sel2 = ' and prop z>'+str(TMD_lower)
            TM = ref.select_atoms(sel1+sel2)
            TM.write(path+names[0]+'.pdb')
        else:
            ref.atoms.write(path+names[0]+'.pdb')
    else:
        for i in range(1, len(names)):
            conf =  path + names[i] 
            top = conf
            u = MDAnalysis.Universe(top, conf,) 
            print('align', names[i], 'with', names[0], 'as reference')
            print('TMD_lower', TMD_lower, 'TMD_higher', TMD_higher)
            align_and_write(universe=u, reference=ref, names=[names[i] ,names[0]],
                            out_path=path, TMD_lower=TMD_lower, TMD_higher=TMD_higher)
    ### hole analysis ###
    aligned_path = path
    fig, ax = plt.subplots()
    plt.title(title, fontsize=f_size)
    colors = ['black', 'blue', 'orange', 'purple','green','red','gray', 'brown',
              'cyan', 'violet', 'olive', 'peru', 'slategray',
    ]
    for i in range(len(names)):
        midpoints, means = hole_analysis(name=names[i]+'.pdb', path=aligned_path, #typ='pdb',
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
    fig.savefig(path + save[:-1] +'HOLE_pathwayprofile.png', bbox_inches='tight')
    #df.to_csv(outpath + '/Pathway_HOLE_comparison_pdb.csv', sep=',', index=False, header=True)
    plt.show()
    ### visualise pathway ###
    pathways = []
    for count, name in enumerate(names):
        #try:
            print(name, '### visualise pathway ###')
            if count<2:
                #! echo $path/{name}.pdb
                #!cp $path/{name}.pdb .
                print('Not copying ...')
            else:
                print('last one..., hard coded')
                path = path + '/PROD/'
                #!cp $path/prod.pdb .
                name = 'prod'
            print('pdbfile=path+name', path+name)
            ha2 =  hole2.hole(
                    pdbfile=path+name ,#+'.pdb',
                    #cpoint='center_of_geometry',
                    executable=hole_exe,
                    #tmpdir=path,
                    #sph_process=sph_proc,
                    sphpdb_file=path+name[:-4]+".sph",
                    end_radius=end_radius,
                    keep_files=True
            )
            #pathway = visualisation(name = names[0], path=path, out=0, end_radius=end_radius )
            #pathways.append(pathway)
            #! mv {name}.sph $path
            #! rm {name}.pdb
        #except:
        #    print('ERROR with', name, 'no SPH file generated')
    #return pathway
    return fig, df
